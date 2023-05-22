import numpy as np
import pandas as pd

def read_b11(filename, index_name = "Z"):
    """
    inputs: 
            filename: str,      name of the file or path to directory
            index_name: str,    name of the index you want to sort on, options: "P", "time", "Z","T","U","TI","PO3","dir","spd"
    outputs:
            data_pandas: pandas DataFrame   dataframe with the data from the b11 file.
    """
    #open the file
    with open(filename) as file:
        binary_data = file.read()

    #split based on this z string that is in every file just before the data starts
    header= binary_data.split("zzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzz")[-1]
    
    #split sentences
    header = header.split("\n")
    
    #make empty data list
    data = []

    #loop through every line with numbers in it
    for i in range(1,len(header)-1):

        #create empty row list
        row = []

        #split based on spaces
        for j in header[i].split(" "):
            
            #check if splits are not empty
            if j !="":

                #append item to row
                row.append(j)

        #append row to list
        data.append(row)
    
    #make the data into a readible pandas file
    data_pandas = pd.DataFrame(data,columns=["P", "time", "Z","T","U","TI","PO3","dir","spd"])
    data_pandas = data_pandas.set_index(index_name)

    #return data
    return data_pandas



print(read_b11("db000106.b11", None))