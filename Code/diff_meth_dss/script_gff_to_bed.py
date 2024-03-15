import argparse
import os

def process(path: str):
    out_path_f = f"{os.path.splitext(path)[0]}_forward.tsv"
    out_path_r = f"{os.path.splitext(path)[0]}_reverse.tsv"
    with open(path, "r") as f, open(out_path_f, "w") as o_f, open(out_path_r, "w") as o_r:
        for _ in range(7):
            next(f)
        for line in f:
            if line.startswith("##"):
                break
            line = line.strip().split("\t")
            chrom, feature, start, end, strand, attribute = line[0], line[2], line[3], line[4], line[6], line[8]
            gene_type = attribute.split(";")[2].strip().replace("gene_type=", "")
            gene_name = attribute.split(";")[3].strip().replace("gene_name=", "")
            gene_id = attribute.split(";")[1].replace("gene_id=", "")
            if strand == "-":
                o_r.write("\t".join([chrom, start, end, strand, feature, gene_type, gene_name, gene_id])+"\n")
            else:
                o_f.write("\t".join([chrom, start, end, strand, feature, gene_type, gene_name, gene_id])+"\n")


if __name__=="__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--path", type=str)

    args = parser.parse_args()
    process(args.path)