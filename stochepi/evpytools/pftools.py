import numpy as np
import scipy.stats as sts
import xml.etree.ElementTree as ET
import warnings

from evpytools import auxiliary as aux


def extract_pfilter_data(xmlfilename, idx=-1, parnames=None):
    """
    Extract data from xml file produced by estavoir (pfilter mode)
    """
    tree = ET.parse(xmlfilename)
    root = tree.getroot()
    
    ## initial parameter values
    initial_params = root.find("initial_params")
    if parnames is None: ## use all unlocked parameters
        if initial_params is None:
            parnames = []
            warnings.warn("element 'initial_params' not present in xml tree")
        else:
            parnames = [x.attrib["name"] for x in initial_params.findall("parameters/param")]
    ## todo: export initial parameter values
            
    
    iterf_steps = root.findall("iterated_filtering_step")

    particle_filters = {
        xs.attrib["id"] : xs.findall("particle_filter")
        for xs in iterf_steps[idx].findall("subject")
    }

    pfIDs = sorted(list(particle_filters.keys()))

    paths = {
        xs.attrib["id"] : xs.findall("path")
        for xs in iterf_steps[idx].findall("subject")
    }

    ranges = {
        k : [ps.findall("pred_range/state") for ps in particle_filters[k]]
        for k in particle_filters.keys()
    }

    pred_medians = {
        k : [ps.find("pred_median/state") for ps in particle_filters[k]]
        for k in particle_filters.keys()
    }

    filter_medians = {
        k : [ps.find("filter_median/state") for ps in particle_filters[k]]
        for k in particle_filters.keys()
    }

    param_ranges = {
        pn : [ps.findall(f"param_range/parameters/param[@name='{pn}']") for ps in iterf_steps]
        for pn in parnames
    }

    param_medians = {
        pn : [ps.find(f"param_median/parameters/param[@name='{pn}']") for ps in iterf_steps]
        for pn in parnames
    }
    
    pf_data_dict = {
        "iterf_steps" : iterf_steps,
        "particle_filters" : particle_filters,
        "pfIDs" : pfIDs,
        "paths" : paths,
        "ranges" : ranges,
        "pred_medians" : pred_medians,
        "filter_medians" : filter_medians,
        "param_ranges" : param_ranges,
        "param_medians" : param_medians,
        "parnames" : parnames,
    }
    
    return pf_data_dict


def makeParDict(r, parnames, medians, ranges):
    def select(par, r):
        if len(par) == 1:
            return float(par[0].attrib["val"])
        else: 
            return float(par[r].attrib["val"])
    ## make the dictionary
    pd = {
        pn : {
            "median" : [select(x, r) for x in medians[pn]],
            "range" : [[select(x, r) for x in xs] for xs in ranges[pn]]
        }
        for pn in parnames
    }       
    return pd
