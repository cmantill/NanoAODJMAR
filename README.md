# NanoAODJMAR

This is a [NanoAOD](https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookNanoAOD) framework for advance developments of jet algorithms. 
The current preliminary content of this development branch can be seen [here](http://algomez.web.cern.ch/algomez/testWeb/JMARnanoplusBTag_content_v01.html) and the size [here](http://algomez.web.cern.ch/algomez/testWeb/JMARnanoplusBTag_size_v01.html).
This format can be used with [fastjet](http://fastjet.fr) directly.

## Recipe

**THIS IS A DEVELOPMENT BRANCH**

For **UL** 2016, 2017 and 2018 data and MC **NanoAODv6** according to the [XPOG](https://gitlab.cern.ch/cms-nanoAOD/nanoaod-doc/-/wikis/Releases/NanoAODv6) and [PPD](https://twiki.cern.ch/twiki/bin/view/CMS/PdmVLegacy2017Analysis) recommendations:

(Using CMSSW_10_6_19 because it has PN already)

```
cmsrel  CMSSW_10_6_19
cd CMSSW_10_6_19/src
cmsenv
git cms-addpkg PhysicsTools/NanoAOD
git cms-addpkg RecoBTag/Combined
git cms-merge-topic andrzejnovak:DDXV2_106
#git cms-merge-topic andrzejnovak:614nosort
git clone https://github.com/cms-data/RecoBTag-Combined.git RecoBTag/Combined/data
git clone https://github.com/cmantill/PFNano.git PhysicsTools/PFNano
#git clone https://github.com/cms-jet/PFNano.git PhysicsTools/PFNano
scram b -j 10
cd PhysicsTools/PFNano/test
```
Note: When running over a new dataset you should check with [the nanoAOD workbook twiki](https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookNanoAOD#Running_on_various_datasets_from) to see if the era modifiers in the CRAB configuration files are correct. The jet correction versions are taken from the global tag.

## Local Usage:

All  python config files were produced with `cmsDriver.py`. 

### How to create python files using cmsDriver

Two imporant parameters that one needs to verify in the central nanoAOD documentation are `--conditions` and `--era`. 
- `--era` options from https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookNanoAOD
- `--conditions` can be found here https://twiki.cern.ch/twiki/bin/view/CMS/PdmV

A an example set of `cmsDriver.py` commands to create configs can be found in `make_configs_preUL.sh`

```
bash make_configs_preUL.sh  # run to produce configs
bash make_configs_preUL.sh  -e # run to actually execute configs on 1000 events
```

## How to create website with nanoAOD content

To create nice websites like [this one](http://algomez.web.cern.ch/algomez/testWeb/JMECustomNano102x_mc_v01.html#Jet) with the content of nanoAOD, use the `inspectNanoFile.py` file from the `PhysicsTools/nanoAOD` package as:
```
python PhysicsTools/NanoAOD/test/inspectNanoFile.py NANOAOD.root -s website_with_collectionsize.html -d website_with_collectiondescription.html
```

## Submission to CRAB

Samples can be submitted to crab using the `submit_all.py` script. Run with `-h` option to see usage. Example can look like this:

```
python submit_all.py -c nano_config.py -s T2_DE_RWTH -f datasets/text_list.txt  -o /store/user/$USER/PFNano/  --ext test --test -d crab_noinpts

```
## Processing data

When processing data, a lumi mask should be applied. The so called golden JSON should be applicable in most cases. Should also be checked here https://twiki.cern.ch/twiki/bin/view/CMS/PdmV

Golden JSON, pre-UL
```
# 2016
jsons/Cert_271036-284044_13TeV_23Sep2016ReReco_Collisions16_JSON.txt
# 2017 
jsons/Cert_294927-306462_13TeV_EOY2017ReReco_Collisions17_JSON_v1.txt
# 2018
jsons/Cert_314472-325175_13TeV_17SeptEarlyReReco2018ABC_PromptEraD_Collisions18_JSON.txt
```

Include in crab submission as `--lumiMask jsons/...txt`

## Documenting the Extended NanoAOD Samples

Please document the input and output datasets on the following twiki: https://twiki.cern.ch/twiki/bin/view/CMS/JetMET/JMARNanoAODv1. For the MC, the number of events can be found by looking up the output dataset in DAS. For the data, you will need to run brilcalc to get the total luminosity of the dataset. See the instructions below. 


## Running brilcalc
These are condensed instructions from the lumi POG TWiki (https://twiki.cern.ch/twiki/bin/view/CMS/TWikiLUM). Also see the brilcalc quickstart guide: https://twiki.cern.ch/twiki/bin/viewauth/CMS/BrilcalcQuickStart.

Note: brilcalc should be run on lxplus. It does not work on the lpc.

Instructions:

1.) Add the following lines to your .bashrc file (or equivalent for your shell). Don't forget to source this file afterwards!

    export PATH=$HOME/.local/bin:/cvmfs/cms-bril.cern.ch/brilconda/bin:$PATH
    export PATH=/afs/cern.ch/cms/lumi/brilconda-1.1.7/bin:$HOME/.local/bin:$PATH
    
2.) Install brilws:

    pip install --install-option="--prefix=$HOME/.local" brilws
    
3.) Get the json file for your output dataset. In the area in which you submitted your jobs:

    crab report -d [your crab directory]
    
The processedLumis.json file will tell you which lumi sections you successfully ran over. The lumi sections for incomplete, failed, or unpublished jobs are listed in notFinishedLumis.json, failedLumis.json, and notPublishedLumis.json. More info can be found at https://twiki.cern.ch/twiki/bin/view/CMSPublic/CRAB3Commands#crab_report.
    
4.) Run brilcalc on lxplus:

    brilcalc lumi -i processedLumis.json -u /fb --normtag /cvmfs/cms-bril.cern.ch/cms-lumi-pog/Normtags/normtag_PHYSICS.json -b "STABLE BEAMS"
    
The luminosity of interest will be listed under "totrecorded(/fb)." You can also run this over the other previously mentioned json files.
    
Note: '-b "STABLE BEAMS"' is optional if you've already run over the golden json. 
        Using the normtag is NOT OPTIONAL, as it defines the final calibrations and detectors that are used for a given run.
