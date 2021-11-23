# GeneplexusPublic
The point of the repo is to a get code base that can be imported by ADS folks without them having to alter any of the python code

The process is split into two parts like the webserver; 1) validate genes and 2) run the model. Example of how to run each is below.

Currently, the repo can be run anywhere in the krishnan lab research space on the HPCC. It also defaults to just printing out some of the inputs that would be used.

### Validate Gene
`python geneplexus.py --task validate --input input_genes.txt`
The input is the file path to the list of input genes (hard coded now, probably need to change, see below)

### Run the model
`python geneplexus.py --task run_model --net_type BioGRID --features Embedding --GSC GO --CV doCV --mg 1000`

Some thoughts and questions

1. Don't know what to do with different input types, right now hard coded in but webserver can handle loading a file and copying in a list
2. How to save outputs or store them for best use in cloud
3. How to get this to run on command line (i.e. geneplexus -t validate) and/or turn into a pip package
4. Can we make a repo/package that is good for cloud and stand alone version. The scripts in `src_hidden` are not used currently but are scripts that can generate all the files at `/mnt/research/compbio/krishnanlab/projects/GenePlexus/repos/GenePlexusBackend/`
5. In the utls.py script, can ignore everythis below `if __name__ == __main__` part, that will probably be removed later
6. Not sure how this compares to the code Doug was actually using. I know the validate part is going be a little different as Doug implemented doing that for all the networks at once himself. I just added that to my script and it was a little bit different.
