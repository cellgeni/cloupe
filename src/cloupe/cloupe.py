import os
import csv
import json
import logging
import struct
import zlib
# import matplotlib.pyplot as plt
# import matplotlib.pyplot as plt
# import matplotlib.image as mpimg


# https://www.10xgenomics.com/datasets/visium-cytassist-gene-and-protein-expression-library-of-human-tonsil-with-add-on-antibodies-h-e-6-5-mm-ffpe-2-standard
# file_path = "datasets/CytAssist_FFPE_Protein_Expression_Human_Tonsil_cloupe.cloupe"

#iget -Kvfr /seq/illumina/runs/49/49384/spaceranger/spaceranger210_count_49384_pSKI_SP15018739_GRCh38-2020-A/cloupe.cloupe
# file_path = "spaceranger210_count_49384_pSKI_SP15018739_GRCh38-2020-A.cloupe"

class Cloupe(object):

    def __init__(self, cloupe_path):
        self.cwd = os.getcwd()
        self.cloupe_path = cloupe_path

        # read main header this is a json that has basic info
        # i'm assuming the header won't be larger than 4096 bytes and it's x00 padded at the end
        self.header = self.read_header()

        # index block will tell us where the data sections start/end
        # it's key for transversing all the sections:
        #   Runs, Matrices, Submatrices, Analyses, Projections, Clusterings, CellTracks,
        #   GeneTracks, Metrics, GeneLists, DiffExps, FeatureClasses, BarcodeClasses,
        #   PeakLocations, FragmentsIndex, GeneAnnotations, TranscriptAnnotations, CellDataTables
        #   BoolFilters, SpatialImages, SpatialImageTiles, SpatialImageSettings, SpatialEnrichments
        #   SpotDeconvolutions, FeatureLinkages, ClonotypeFiles
        self.index_block = self.read_block(
            start=self.header["indexBlock"]["Start"],
            end=self.header["indexBlock"]["End"],
            as_json=True,
        )

        # read Runs
        self.runs = []
        for run in self.index_block.get("Runs", []):
            if "Metadata" in run:
                run["Metadata"] = self.read_block(
                    start=run["Metadata"]["Start"],
                    end=run["Metadata"]["End"],
                    as_json=True,
                )
            self.runs.append(run)

        # read Metrics
        self.metrics = []
        for metric in self.index_block.get("Metrics", []):
            self.metrics.append(
                self.read_block(
                    start=metric["Contents"]["Start"],
                    end=metric["Contents"]["End"],
                    as_json=True,
                )
            )

        # read Analyses
        self.analyses = []
        for analysis in self.index_block["Analyses"]:
            self.analyses.append(analysis)

        # read Projections
        self.projections = {}
        for projection in self.index_block["Projections"]:
            try:
                # read data block
                pblock_bytes = self.read_block(
                    start=projection["Matrix"]["Start"],
                    end=projection["Matrix"]["End"],
                )
                pblock_count = projection["Matrix"]["ArraySize"]
                pblock = struct.unpack(f"{pblock_count}d", pblock_bytes)

                # read index block
                pindex_bytes = self.read_block(
                    start=projection["Matrix"]["Index"]["Start"],
                    end=projection["Matrix"]["Index"]["End"],
                )
                pindex_count = projection["Matrix"]["Index"]["ArraySize"]
                pindex = struct.unpack(f"{pindex_count}Q", pindex_bytes)

                # reashape so it has the dimensions of the matrix we want
                cols, rows = tuple(projection["Dims"])
                self.projections[projection["Name"]] = [
                    pblock[x : x + rows] for x in range(0, len(pblock), rows)
                ]
            except Exception as ex:
                logging.error(f"Projection '{projection['Name']}' error: {ex}")

        # read Matrices
        # Matrices have lots of properties, i'm only picking some, but the complete list is:
        #   Name', Uuid, FormatVersion, ParentUuid, Reference, GeneCount, CellCount, Genes, GeneNames,
        #   FeatureCount, BarcodeCount, Metadata, Barcodes, FeatureIds, FeatureNames, FeatureSecondaryNames,
        #   FeatureTypeOffsets, FeatureTypeIndices, CSCValues, CSCPointers, CSCIndices, CSRValues, CSRPointers,
        #   CSRIndices, UMICounts
        self.matrices = []
        for matrix in self.index_block["Matrices"]:
            try:
                for prop in ["Barcodes","FeatureIds","FeatureNames","FeatureSecondaryNames",]:
                    if prop not in matrix:
                        continue
                    mblock = self.read_block(
                        start=matrix[prop]["Start"], end=matrix[prop]["End"]
                    )
                    step = matrix[prop]["ArrayWidth"]
                    stop = matrix[prop]["ArraySize"] * step
                    matrix[prop] = [
                        mblock[i : i + step]
                        .decode("utf-8", errors="replace")
                        .strip("\x00")
                        for i in range(0, stop, step)
                    ]
                umicounts_bytes = self.read_block(
                    start=matrix["UMICounts"]["Start"],
                    end=matrix["UMICounts"]["End"],
                )
                umicounts_count = matrix["UMICounts"]["ArraySize"]
                matrix["UMICounts"] = struct.unpack(
                    f"{umicounts_count}Q", umicounts_bytes
                )
                self.matrices.append(matrix)
            except Exception as ex:
                logging.error(f"Matrix '{matrix['Name']}' error: {ex}")

        # read the next header
        # this looks like that is where the custom user data is written
        # there could be several chained headers by nextHeaderOffest but i'm only
        # gonna parse the first one and i'll only get CellTracks, however, more
        # stuff could be here that we want to export
        self.celltracks = []
        if self.header["nextHeaderOffset"] != self.get_file_size():
            self.next_header = self.read_header(self.header["nextHeaderOffset"])
            self.next_index_block = self.read_block(
                start=self.next_header["indexBlock"]["Start"],
                end=self.next_header["indexBlock"]["End"],
                as_json=True,
            )
            # now that we have the index block we get the 'userCreated' celltracks
            self.celltracks = []
            for celltrack in self.next_index_block["CellTracks"]:
                celltrack["Metadata"] = self.read_block(
                    start=celltrack["Metadata"]["Start"],
                    end=celltrack["Metadata"]["End"],
                    as_json=True,
                )
                vblock_bytes = self.read_block(
                    start=celltrack["Values"]["Start"],
                    end=celltrack["Values"]["End"],
                )
                vblock_count = celltrack["Values"]["ArraySize"]
                groups = celltrack["Metadata"]["groups"]
                group_keys = struct.unpack(f"{vblock_count}h", vblock_bytes)
                celltrack["Values"] = [groups[gk] if gk!=-1 else "UNASSIGNED" for gk in group_keys]
                self.celltracks.append(celltrack)

    def read_header(self, position=0):
        with open(self.cloupe_path, "rb") as f:
            f.seek(position)
            header = f.read(4096)
            return json.loads(header.decode("utf-8", errors="replace").strip("\x00"))

    def get_file_size(self):
        with open(self.cloupe_path, "rb") as f:
            f.seek(0, 2)
            size = f.tell()
            return size

    def __repr__(self):
        return str(self.index_block["Runs"])

    # trick for partial object decompression
    def decompress(self, chunk):
        gz = zlib.decompressobj(31)
        return gz.decompress(chunk)

    # read cloupe block
    def read_block(self, start, end, as_json=False):
        with open(self.cloupe_path, "rb") as f:
            f.seek(start)
            read_size = end - start
            block = f.read(read_size)
            if block[:2] == b"\x1f\x8b":
                logging.debug("compressed block")
                block = self.decompress(block)
            if as_json:
                return json.loads(block.decode("utf-8", errors="replace"))
            return block

    # writing the barcodes to the csv file
    def barcodes_writer(self):
        barcodes = self.cwd + '/barcodes.csv'
        with open(barcodes, 'w', newline="") as barcodes_file:
            # Initialize the CSV writer
            csv_writer = csv.writer(barcodes_file)
            csv_writer.writerow(["Barcodes"])
            csv_writer.writerows(
                    [[element] for element in self.matrices[0]['Barcodes']]
            )  # Efficient one line code

    # write the features to the csv file
    def features_writer(self):
        features = self.cwd + '/features.csv'
        with open(features, 'w', newline="") as features_file:
            csv_writer = csv.writer(features_file)
            csv_writer.writerow(["FeatureIds", "FeatureNames"])
            csv_writer.writerows([[str(element) for element in pair]  # Remove parentheses, single quotes, and spaces
                        for pair in zip(
                        self.matrices[0]["FeatureIds"],
                        self.matrices[0]["FeatureNames"],
                    )])# Using zip command along with list comprehension to convert the two list into strings for the csv as elements

    # write annotations
    def annotations_writer(self):
        annotations = self.cwd + '/annotations.csv'
        with open(annotations, 'w', newline="") as annotations_file:
            csv_writer = csv.writer(annotations_file)
            csv_writer.writerow(["Barcodes"]+[element["Name"] for element in self.celltracks])
            csv_writer.writerows(
                [
                    [barcode, *values]
                    for barcode, *values in zip(
                    self.matrices[0]["Barcodes"],
                    *(track["Values"] for track in self.celltracks)
                )
                ]
            )

    # plotting the spatial info
    def spatial_projection(self):
        spatial_plot = self.cwd + '/projection_spatial.csv'
        spatial_embedding = self.projections['Spatial']
        projections_1 = [["Barcodes"]+self.matrices[0]["Barcodes"]]
        for column in range(len(spatial_embedding)):
            projections_1.append(["d{}".format(column)]+[element for element in spatial_embedding[column]])
        transposed_projections_1 = zip(*projections_1)
        with open(spatial_plot, "w", newline="") as spatial_projection_file:
            csv_writer = csv.writer(spatial_projection_file)
            csv_writer.writerows(transposed_projections_1)

        # plt.scatter(x=spatial_embedding[0], y=spatial_embedding[1])
        # plt.gca().invert_yaxis()
        # plt.show()

    # plotting the tsne info
    def tsne_projection(self):
        tsne_plot = self.cwd + '/projection_tsne.csv'
        tsne_embedding = self.projections['tsne']
        projections_2 = [["Barcodes"] + self.matrices[0]["Barcodes"]]
        for column in range(len(tsne_embedding)):
            projections_2.append(["d{}".format(column)] + [element for element in tsne_embedding[column]])
        transposed_projections_2 = zip(*projections_2)
        with open(tsne_plot, "w", newline="") as tsne_projection_file:
            csv_writer = csv.writer(tsne_projection_file)
            csv_writer.writerows(transposed_projections_2)


if __name__=="__main__":
    logging.basicConfig(level=logging.INFO, format="[%(asctime)s][%(levelname)s] %(message)s")
    cloupe_object = Cloupe("/Users/nj9/Downloads/spaceranger210_count_49384_pSKI_SP15018739_GRCh38-2020-A.cloupe")
    # cloupe_object.barcodes_writer()
    # cloupe_object.features_writer()
    # cloupe_object.annotations_writer()
    # cloupe_object.spatial_projection()
# Bye-bye



