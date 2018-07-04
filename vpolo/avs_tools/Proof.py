# import click
from EqClass import Eqclasses
import Alevin
import Utools
import Kallisto
import CellRanger
import Alevin_MST
import sys

def plot_predictions(alv, utl, kal, crn, alv_mst):
    print("\n\n\n\n===================")
    print("IGNORE logs above this\n===================")
    print("Alevin Predictions: ")
    print(alv)
    print("Alevin_MST Predictions: ")
    print(alv_mst)
    print("Utools Predictions: ")
    print(utl)
    print("Kallisto Predictions: ")
    print(kal)
    print("CellRanger Predictions: ")
    print(crn)

def run(eq_file):
    '''
    Run tests for all the methods
    :param eq_file: path of the file having eqclass structure
    :return: None
    '''
    eq_class_obj = Eqclasses(eq_file)

    # get gene level prediction dictionary for each tool
    utl_prediction = Utools.get_prediction(eq_class_obj)
    alv_prediction = Alevin.get_prediction(eq_class_obj)
    alv_mst_prediction = Alevin_MST.get_prediction(eq_class_obj)
    kal_prediction = Kallisto.get_prediction(eq_class_obj)
    crn_prediction = CellRanger.get_prediction(eq_class_obj)

    # plot the histogram of the prediction for each tool
    plot_predictions(alv_prediction, utl_prediction,
                     kal_prediction, crn_prediction,
                     alv_mst_prediction)

if __name__ == "__main__":
    run(sys.argv[1])
