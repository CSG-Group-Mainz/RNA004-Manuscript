import argparse
import os
def prep_pileup(path: str):
    out_path_forw = f"{os.path.splitext(path)[0]}_forw_prepped.bed"
    out_path_rev = f"{os.path.splitext(path)[0]}_rev_prepped.bed"

    with open(path, "r") as f, open(out_path_forw, "w") as o_f, open(out_path_rev, "w") as o_r:
        o_f.write("chr\tpos\tN\tX\n")
        o_r.write("chr\tpos\tN\tX\n")
        for line in f:
            line = line.strip().split("\t")
            line[9] = line[9].split(" ")

            if int(line[9][0]) > 4:
                new_line = "\t".join([line[0], line[2], line[9][0], line[9][2]])
                if line[5] == "-":
                    o_r.write(new_line+"\n")
                else:
                    o_f.write(new_line+"\n")

if __name__=="__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--path", type=str)

    args = parser.parse_args()
    prep_pileup(args.path)
