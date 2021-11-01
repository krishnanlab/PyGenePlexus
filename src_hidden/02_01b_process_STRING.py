import numpy as np
import pickle
import os
import time


def load_convert_dict():
    fp_dict = '/mnt/research/compbio/krishnanlab/data/MyGeneInfo/20201029_Entrez_Multiple-Species/ID_conversions/'
    with open(fp_dict+'Homo_sapiens__ENSP-to-Entrez__All-Mappings.pickle', 'rb') as mydict:
        convert_dict = pickle.load(mydict)
    return convert_dict
    
def check_exists(FileNameList):
    for item in FileNameList:
        if os.path.isfile(item) == True:
            raise FileExistsError('{} already exists'.format(item))
            
def main(File_org):
    Com_dict = {}
    Exp_dict = {}
    with open(File_org, mode='r') as f_old:
        for idx, line in enumerate(f_old):
            if idx == 0:
                continue
            else:
                # get the uniprot IDs and final and initial scores
                ENSPa = line.split(' ')[0].split('.')[1]
                ENSPb = line.split(' ')[1].split('.')[1]
                ExpW = int(line.split(' ')[6])
                ComW = int(line.split(' ')[9])
            # convert ENSP to entrez
            try:
                EntAs = convert_dict[ENSPa]
                EntBs = convert_dict[ENSPb]
                for EntA in EntAs:
                    for EntB in EntBs:
                        EntA = int(EntA)
                        EntB = int(EntB)
                        # sort EntA and EntB and remove self edges
                        if EntA == EntB:
                            continue
                        elif EntB < EntA:
                            EntB, EntA = EntA, EntB
                        else:
                            pass
                        mykey = '{}_{}'.format(EntA,EntB)
                        if ComW > 0:
                            if mykey not in Com_dict:
                                Com_dict[mykey] = ComW
                            else:
                                prev_value = Com_dict[mykey]
                                if ComW > prev_value:
                                    Com_dict[mykey] = ComW
                                else:
                                    pass
                        if ExpW > 0:
                            if mykey not in Exp_dict:
                                Exp_dict[mykey] = ExpW
                            else:
                                prev_value = Exp_dict[mykey]
                                if ExpW > prev_value:
                                    Exp_dict[mykey] = ExpW
                                else:
                                    pass
            except KeyError:
                continue

    return Com_dict, Exp_dict
                        
def make_edgelist(Com_dict,Exp_dict,File_Com,File_Exp): 
    for item in [(Com_dict,File_Com),(Exp_dict,File_Exp)]:
        mylist = [(int(key.split('_')[0]),int(key.split('_')[1]),float(val)) for key, val in item[0].items()]
        mylist = sorted(mylist, key=lambda element: (element[0], element[1]))
        myarr = np.array(mylist)
        myarr[:,2] = myarr[:,2]/1000.
        np.savetxt(item[1],myarr,fmt=['%i','%i','%.3f'],delimiter='\t')


if __name__ == "__main__":
    
    # This will make the edgelist from anything that has a combined score
    # or anything that was experimental obtained in just human experiments
    # make paths to the files
    fp = '/mnt/research/compbio/krishnanlab/projects/GenePlexus/repos/GenePlexusBackend/data_hidden2/Network_Downloads/STRING/'
    File_org = fp + '9606.protein.links.detailed.v11.0.txt'
    fp2 = '/mnt/research/compbio/krishnanlab/projects/GenePlexus/repos/GenePlexusBackend/data_backend2/Edgelists/'
    File_Com = fp2 + 'STRING.edg'
    File_Exp = fp2 + 'STRING-EXP.edg'
    # check if the modified files exist
    # FileList = [File_Com,File_Exp]
    # check_exists(FileList)
    #load the conversion dictionary
    convert_dict = load_convert_dict()
    # # run main file
    print('Making edgelist dictionary')
    tic = time.time()
    Com_dict, Exp_dict = main(File_org)
    print('The numner of seconds it took to run last step is',int(time.time()-tic))
    # # make edgelist
    print('Making final edgelist')
    tic = time.time()
    make_edgelist(Com_dict,Exp_dict,File_Com,File_Exp)
    print('The numner of seconds it took to run last step is',int(time.time()-tic))

    
