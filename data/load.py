# -*- coding: utf-8 -*-
import pathlib, json

def load_json(filename):
    
    path = pathlib.Path(__file__).parent / (filename + ".json")
        
    with open(path) as file_json:
        file_dict = json.load(file_json)
        
    return file_dict


def data(filename, verbose=False):
    
    if verbose:
        info(filename)
    
    file_dict = load_json(filename)
    layers = None
    
    for key in list(file_dict.keys()):
        val = file_dict[key]
        if isinstance(val, str):
            del file_dict[key]
    
    for key, val in file_dict.items():
        if isinstance(val, list):
            file_dict[key] = val[0] * val[1]
            
    if "layers" in file_dict:
        for layer in file_dict["layers"].values():
            for param, val in layer.items():
                if isinstance(val, list):
                    layer[param] = val[0] * val[1]
        layers = list(file_dict["layers"].values())
        del file_dict["layers"]
                
    return file_dict, layers


def info(filename):
    file_dict = load_json(filename)
    
    max_title_length = max([len(key)+len(str(val)) for key, val in file_dict.items() 
                            if not isinstance(val, dict)])
    max_info_length = max([len(key) for key, val in file_dict.items() 
                           if not isinstance(val, dict)])
    max_data_length = max([len(key) for key, val in file_dict.items() 
                           if not isinstance(val, str)])
    
    print(f"\n--- MACHINE INFORMATIONS ---{(max_title_length-23) * '-'}")
    for key, val in file_dict.items():
        if isinstance(val, str):
            print(f"{key}: {(max_info_length - len(str(key))) * ' '} {val}")
            if key == "Type": print("")
    
    print(f"\n--- MACHINE DATA ---{(max_title_length-15) * '-'}")           
    for key, val in file_dict.items():
        if isinstance(val, int):
            print(f"{key}{(max_data_length - len(str(key))) * ' '} = {val}")
    
    for key, val in file_dict.items():
        if isinstance(val, list):
            print(f"{key}{(max_data_length-len(str(key)))*' '} = {val[0]}{(5-len(str(val[0])))*' '}{val[-1]}")
            
    if "layers" in file_dict:
        print(f"\n--- LAYER DATA ---{(max_title_length-13) * '-'}")
        for lay_no, lay_info in file_dict["layers"].items():
            print(f"--- Layer: {lay_no} ---")
            max_key_length = max([len(key) for key in lay_info])
            for key, val in lay_info.items():
                if isinstance(val, list):
                    print(f"{key}{(max_key_length  - len(str(key))) * ' '} = {val[0]} {val[-1]}")
                elif isinstance(val, int):
                    print(f"{key}{(max_key_length  - len(str(key))) * ' '} = {val}")
                else:
                    print(f"{key}{(max_key_length  - len(str(key))) * ' '}:  {val}")
            print("")
    print(f"----------------------------{(max_title_length-23) * '-'}\n")    
    
            
            
            
            
