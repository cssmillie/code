library(argparse)

if(interactive()){

    # set default arguments
    # ---------------------
    args = list()

} else {

    # read input arguments
    # --------------------
    parser = ArgumentParser()
    parser$add_argument('-i', help='input file')
    parser$add_argument('-o', help='output file')
    args = parser$parse_args()

}
