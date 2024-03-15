import argparse
import os

def process(path: str):
    out_path = f"{os.path.splitext(path)[0]}.bed"
    with open(path, "r") as f, open(out_path, "w") as o:
        next(f)
        for line in f:
            line = line.strip().split("\t")
            chrom, end  = line[0], int(line[1])
            o.write("\t".join([chrom, str(end-1), str(end)]+line[2:])+"\n")

if __name__=="__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--path", type=str)

    args = parser.parse_args()
    process(args.path)
