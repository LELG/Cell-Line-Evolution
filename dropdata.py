from __future__ import print_function
import csv


def drop(subpop, results_dir, label):
    ## print info for multiple populations
    ## print by colour of population
    #summary for all runs - add value and overwrite

    cols_sizes = subpop.tree_to_list("size_by_col")

    cols_sizes.sort()

    colours = ['LightGray', 'Gray', 'LightSkyBlue', 'RoyalBlue',
               'LightSeaGreen', 'MediumSeaGreen', 'Khaki', 'Goldenrod',
               'Tomato', 'PaleVioletRed', 'DarkViolet']

    col_totals = [0 for col in colours]

    #merge colours, super lazynaive method
    #will eventually change whole color sys to be dynamic

    for col, size in cols_sizes:
        col_totals[colours.index(col)] += size

    drop_fpath = "{0}/{1}dropdata.csv".format(results_dir, label)
    drop_file = open(drop_fpath, 'a')
    drop_writer = csv.writer(drop_file)
    drop_writer.writerow(colours)
    drop_writer.writerow(col_totals)
    drop_file.close()


def read_drop(test_dir):
    file_mid = open(test_dir+"/sum_hetpop_middist.csv", 'a')
    file_end = open(test_dir+"/sum_hetpop_enddist.csv", 'a')
    #MID
    print("reading...")
    colours = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
    #Read the info for all colour populations, and graph
    filename = "{0}/middropdata.csv".format(test_dir)
    mid_infile = open(filename)
    mid_inreader = csv.reader(mid_infile)
    #sum total counts
    _header = next(mid_inreader)
    for line in mid_inreader:
        sizes = [int(z) for z in line]
        print("and then")
        colours = [x + y for x, y in zip(colours, sizes)]
    print("ALL MID COLOURS ", colours, file=file_mid)
    mid_infile.close()
    #then plot it

    #END
    print("reading...")
    colours = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
    filename = "{0}/enddropdata.csv".format(test_dir)
    end_infile = open(filename)
    end_inreader = csv.reader(end_infile)
    #sum total counts
    _header = next(end_inreader)
    for line in end_inreader:
        sizes = [int(z) for z in line]
        colours = [x + y for x, y in zip(colours, sizes)]
    print("ALL END COLOURS ", colours, file=file_end)
    #then plot it

#parser = argparse.ArgumentParser()
#parser.add_argument('-f','--filez', default='filez')
#opt = parser.parse_args()
#read_drop(opt.filename,'mid')
#read_drop(opt.filename,'end')
