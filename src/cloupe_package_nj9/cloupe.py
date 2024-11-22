import os
import re
import json
import logging
import struct
import zlib
import argparse
import matplotlib.pyplot as plt
import matplotlib.pyplot as plt
import matplotlib.image as mpimg


# https://www.10xgenomics.com/datasets/visium-cytassist-gene-and-protein-expression-library-of-human-tonsil-with-add-on-antibodies-h-e-6-5-mm-ffpe-2-standard
# file_path = "datasets/CytAssist_FFPE_Protein_Expression_Human_Tonsil_cloupe.cloupe"

#iget -Kvfr /seq/illumina/runs/49/49384/spaceranger/spaceranger210_count_49384_pSKI_SP15018739_GRCh38-2020-A/cloupe.cloupe
# file_path = "spaceranger210_count_49384_pSKI_SP15018739_GRCh38-2020-A.cloupe"

class Cloupe(object):

    def __init__(self, cloupe_path):
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


if __name__=="__main__":
    logging.basicConfig(level=logging.INFO, format="[%(asctime)s][%(levelname)s] %(message)s")
    # Create the parser
    parser = argparse.ArgumentParser(description='Cloupe parser which extracts Barcodes, Features and Annotations info')
    # Add arguments
    parser.add_argument('input_file', type=str, help='Path to the input file')
    # Parse the arguments
    args = parser.parse_args()
    # Access the arguments
    # Load cloupe file
    cwd = os.getcwd()
    # hardcoding the location of the files that are generated out of the script
    barcodes = cwd + '/barcodes.csv'
    features = cwd + '/features.csv'
    barcode_clusters = cwd + '/annotations.csv'
    cloupe = Cloupe(args.input_file)

    count_values = {}
    # Getting the count of Features and Barcodes to set the limits for the columns in the csv files
    for prop in ['FeatureCount', 'BarcodeCount']:
        count_values[prop] = cloupe.matrices[0][prop]


    collection = []
    # writing the barcodes to the csv files
    with open(barcodes, 'w') as barcodes_file:
        barcodes_file.write("Barcodes\n")
        barcodes_file.write(
            "\n".join(
            cloupe.matrices[0]['Barcodes'][:count_values['BarcodeCount']]
            )
        )# Efficient one line code

    with open(features, 'w') as features_file:
        features_file.write("FeatureIds,FeatureNames\n")
        features_file.write(
            "\n".join(
                [
                    re.sub(r"[()' ]", "", str(pair))  # Remove parentheses, single quotes, and spaces
                    for pair in zip(
                    cloupe.matrices[0]["FeatureIds"][: count_values["FeatureCount"]],
                    cloupe.matrices[0]["FeatureNames"][: count_values["FeatureCount"]],
                )
                ]
            )
        )# Using zip command along with list comprehension to convert the two list into strings for the csv as elements

    # Plotting the info
    spatial_embedding = cloupe.projections['Spatial']
    plt.scatter(x=spatial_embedding[0], y=spatial_embedding[1])
    plt.gca().invert_yaxis()
    plt.show()

    # Final csv addition
    with open(barcode_clusters, 'w') as barcode_clusters_file:
        barcode_clusters_file.write("Barcodes,{},{}\n".format(cloupe.celltracks[0]["Name"],cloupe.celltracks[1]["Name"]))
        barcode_clusters_file.write(
            "\n".join(
                [
                    re.sub(r"[()' ]", "", str(pair))  # Remove parentheses, single quotes, and spaces
                    for pair in zip(
                    cloupe.matrices[0]["Barcodes"],
                    cloupe.celltracks[0]["Values"],
                    cloupe.celltracks[1]["Values"],
                )
                ]
            )
        )

# Bye-bye



