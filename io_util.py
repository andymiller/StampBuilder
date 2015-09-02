"""
Some IO functions and some wrappers to ease laptop/nersc use
"""
import sys, os
import numpy as np

# determine platform from paths (HACK CITY!)
if os.path.exists("/project/projectdirs/cosmo/data/sdss"):
    print "ON NERSC, using local SDSS data"
    PLATFORM   = "NERSC"

    # set up environment variables for DR9 object
    os.environ["PHOTO_RESOLVE"] = \
        "/project/projectdirs/cosmo/data/sdss/pre13/eboss/resolve/2013-07-29"
    os.environ["BOSS_PHOTOOBJ"] = \
        "/project/projectdirs/cosmo/data/sdss/pre13/eboss/photoObj.v5b"
    os.environ["PHOTO_REDUX"] = ""
else:
    print "Local - running sdss out of scratch/blobs/"
    PLATFORM = 'local'
    sys.path.insert(1, "lib/astrometry.net")
    sys.path.insert(1, "lib/tractor")


# safe import of tractor now
import tractor.sdss as tsdss


#def platform_agnostic_sdss():
#    """ Returns a Tractor SDSS object for handling file io """
#    if PLATFORM == "NERSC":
#        sdss = tsdss.DR9(basedir="/project/projectdirs/das/acmiller/sdss_base/")
#        sdss.useLocalTree()
#        return sdss
#    elif PLATFORM == "local":
#        local_sdss_basedir = 'scratch/blobs'
#        sdss = tsdss.DR9(basedir = local_sdss_basedir)
#        sdss.saveUnzippedFiles(local_sdss_basedir)
#        return sdss
#    else:
#        raise Exception("Platform %s has no sdss directory!"%PLATFORM)
#
def catalog_sdss():
    """ returns a tractor sdss that is pointed at catalog files """
    if PLATFORM == "NERSC":
        sdss = tsdss.DR9(basedir="/project/projectdirs/das/acmiller/sdss_base/")
        sdss.useLocalTree(photoObjs = "/project/projectdirs/cosmo/data/sdss/pre13/eboss/photoObj.v5b",
                          resolve   = "/project/projectdirs/cosmo/data/sdss/pre13/eboss/resolve/2013-07-29")
        return sdss
    elif PLATFORM == "local":
        local_sdss_basedir = 'scratch/blobs'
        sdss = tsdss.DR9(basedir = local_sdss_basedir)
        sdss.saveUnzippedFiles(local_sdss_basedir)
        return sdss
    else:
        raise Exception("Platform %s has no sdss directory!"%PLATFORM)


def photo_sdss():
    if PLATFORM == "NERSC":
        sdss = tsdss.DR9(basedir="/global/projecta/projectdirs/sdss/data/sdss/dr9/boss/photo/redux") #basedir="/project/projectdirs/das/acmiller/sdss_base/")
        sdss.useLocalTree(
            photoObjs = "/global/projecta/projectdirs/sdss/data/sdss/dr9/boss/photoObj",
            resolve   = "/global/projecta/projectdirs/sdss/data/sdss/dr9/boss/resolve/2010-05-23")
        os.environ["PHOTO_REDUX"] = ""
        return sdss
    elif PLATFORM == "local":
        local_sdss_basedir = 'scratch/blobs'
        sdss = tsdss.DR9(basedir = local_sdss_basedir)
        sdss.saveUnzippedFiles(local_sdss_basedir)
        return sdss
    else:
        raise Exception("Platform %s has no sdss directory!"%PLATFORM)


def default_output_dir():
    """ on nersc, store in /project/projectdirs/das/acmiller"""
    if PLATFORM == "NERSC":
        return "/project/projectdirs/das/acmiller"
    elif PLATFORM == "local":
        return "stamps"
    else:
        raise Exception("Platform %s not supported!"%PLATFORM)


def load_mcmc_chains(chain_file_0, num_chains=4, burnin=2500):
    """ Load galaxy mcmc chains from chain files, aggregate """
    src_samp_chains = []
    ll_samp_chains  = []
    eps_samp_chains = []
    for i in range(num_chains):
        chain_i_file = chain_file_0.replace("chain_0", "chain_%d"%i)
        if not os.path.exists(chain_i_file):
            print "chain_file: %s does not exist"%chain_file_0
            continue
        samp_dict = load_samples(chain_i_file)
        if samp_dict['srcs'][burnin] == samp_dict['srcs'][burnin+1]:
            print "chain_file: %s has zeros?"%chain_file_0
            print samp_dict['srcs'][burnin:(burnin+2)]
            continue

        src_samp_chains.append(samp_dict['srcs'][burnin:,0])
        ll_samp_chains.append(samp_dict['ll'][burnin:])
        eps_samp_chains.append(samp_dict['epsilon'][burnin:])
    return src_samp_chains, ll_samp_chains, eps_samp_chains


# write samps to file
SAMPLE_FIELDS = ['epsilon', 'srcs', 'll']
def save_samples(samp_dict, fname):
    f = file(fname, "wb")
    for s in SAMPLE_FIELDS:
        np.save(f, samp_dict[s])
    f.close()

def load_samples(fname):
    f = file(fname, "rb")
    samp_dict = {}
    for s in SAMPLE_FIELDS:
        samp_dict[s] = np.load(f)
    return samp_dict

