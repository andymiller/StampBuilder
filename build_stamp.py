"""
assembles dataset of galaxy stamps from an object file.

The object file has information from a list of sources.  An example SQL
QUERY to build a fits file of object information (submit to CasJob tool):
    select
        *
    into mydb.s82
    from
      dr12.PhotoPrimary
    where
        ra  between 5 and 6 and
        dec between 0 and 1 and
        (run = 106 or run = 206)
    -> s82.fits

Ref:
Dustin's patch script:
  This includes the query to the SQL database for Galaxy Information 
  https://github.com/dstndstn/tractor/blob/master/projects/inference/testblob2.py

Galaxy Information database:
  Where to get galaxy source information (RA, DEC) to build patches
  - http://skyserver.sdss.org/casjobs/login.aspx
  - schema info: http://skyserver.sdss.org/dr12/en/help/browser/browser.aspx#&&history=description+Galaxy+V

"""
import matplotlib
matplotlib.use('Agg')
import pylab as plt
import numpy as np
import sys, os
import fitsio
import io_util
import astrometry.util.fits         as aufits
import astrometry.util.util         as autil
import astrometry.util.plotutils    as aplot
import astrometry.util.resample     as aresample
import astrometry.libkd.spherematch as asphere
import astrometry.sdss              as asdss
import astrometry.util.multiproc    as amultiproc
import tractor
import tractor.sdss                 as tsdss


#########################################################################
# set up arguments for CLI
#########################################################################
import argparse
parser = argparse.ArgumentParser(description='Build a dataset of postage stamps')
parser.add_argument('--output_dir', 
    action  = "store", 
    dest    = "output_dir", 
    default = os.path.join(io_util.default_output_dir(), "stamps"))
parser.add_argument("--source_file", action  = "store",
                                     dest    = "source_file",
                                     default = "fits_files/s82-1k.fits")
parser.add_argument("--idx_keep", action  = "store",
                                  dest    = "idx_keep",
                                  default = 1)
parser.add_argument("--num_proc", action = "store",
                                  dest   = "num_proc",
                                  default = 8)
parser.add_argument('--version', action='version', version='0.1')
results = parser.parse_args()


#######################################################################
# GLOBAL PARAMS i shouldn't have to pass around
#######################################################################
pixscale  = 0.396
pixradius = 25
bands     = 'ugriz'
radius    = pixradius * pixscale / 3600.
srcband   = 'r'
Lanczos   = 3
source_idx = int(results.idx_keep)

def main():
    # load in parsed args
    object_file = results.source_file
    outdir      = results.output_dir
    idx_keep    = int(results.idx_keep)
    Nproc       = results.num_proc
    print """
    =================================================================
    BUILD DATASET CALLED:
      object_file = {object_file}
      outdir      = {outdir}
      idx_keep    = {idx_keep}
      Nproc       = {Nproc}
    =================================================================
     """.format(object_file = object_file,
                outdir = outdir,
                idx_keep = idx_keep,
                Nproc = Nproc)

    # load and clean objects
    T = load_and_clean_objects(object_file, outdir=outdir, idx_keep=idx_keep)

    # write object cat files - will correspond to each stamp
    #write_cat_files(T, outdir=outdir)
    cutToPrimary = False

    # grab the necessary information for each source
    stars = [(ra,dec,[],cutToPrimary,outdir) for ra,dec in zip(T.ra, T.dec)]

    # Write out stamps and cat files
    mp = amultiproc.multiproc(1)
    mp.map(_bounce_one_blob, stars)


def _bounce_one_blob(args):
    try:
        oneblob(*args)
    except:
        print 'Error running oneblob:'
        import traceback
        traceback.print_exc()
        print


def oneblob(ra, dec, addToHeader, cutToPrimary, outdir):
    """ Given an RA, DEC, and outdir ...
         - cuts out a small patch around the source
         - resamples all test blobs to a common pixel grid
         - saves the resulting stamps
    """
    # Identify stamp-ification using the r-band png
    plotfn = os.path.join(outdir, 'stamps-%.4f-%.4f.png' % (ra, dec))
    if os.path.exists(plotfn):
        print '\n======================================================'
        print 'Exists:', plotfn
        print '========================================================\n'
        return []

    # compute stamp output information - (e.g. radius around the given RA, DEC)
    W,H = pixradius*2+1, pixradius*2+1
    targetwcs = autil.Tan(ra, dec, pixradius+1, pixradius+1,
                          -pixscale/3600., 0., 0., pixscale/3600., W, H)

    # get the fields that are in this range
    print """
    ===================================================
    Getting overlapping Run, Camcol, Field values for
    ra, dec = {ra}, {dec}
    """.format(ra=ra, dec=dec)
    RCF = get_overlapping_run_camcol_field_rerun_301(ra, dec,
      io_util.catalog_sdss())

    # create source info table fields found above, write to catfn file
    print """
    ===================================================
    Creating source table from fields for
    ra, dec = {ra}, {dec}
    """.format(ra=ra, dec=dec)
    T = create_source_table_from_fields(RCF, ra, dec,
                                        cutToPrimary, srcband, 
                                        io_util.catalog_sdss())
    catfn = os.path.join(outdir, 'cat-%.4f-%.4f.fits' % (ra,dec))
    T.writeto(catfn)

    # track output files (i think this will just be [catfn])
    outfns = []
    outfns.append(catfn)

    # construct a multi-stamp image for each band
    for band in bands:

        # For each band, write out a fits image with a fit PSF in the header
        print """

            ===================================================
            Making resampled PSF Images for
                ra, dec = {ra}, {dec}
                band    = {band}
        """.format(ra=ra, dec=dec, band=band)
        print "photo redux environ?", os.environ["PHOTO_REDUX"]
        resampled_imgs = make_resampled_psf_images(
            RCF, band, ra, dec, io_util.photo_sdss(), targetwcs, W, H, 
            addToHeader)

        # write out a single FITS file 
        fn = stamp_filename(outdir, band, ra, dec)
        print 'writing', fn
        clobber = True
        for img, iv, hdr in resampled_imgs:
            fitsio.write(fn, img.astype(np.float32), clobber=clobber, header=hdr)
            fitsio.write(fn, iv.astype(np.float32))
            if clobber:
                outfns.append(fn)
            clobber = False

        # create stamps image for the one
        if band == 'r':
            plt.figure(figsize=(8,8))
            plt.subplots_adjust(left=0.01, right=0.99, bottom=0.01, top=0.99,
                                hspace=0.05, wspace=0.05)
            N     = len(resampled_imgs)
            ncols = int(np.ceil(np.sqrt(float(N))))
            nrows = int(np.ceil(float(N) / ncols))
            plt.clf()
            for k, (img, iv, hdr) in enumerate(resampled_imgs):
                plt.subplot(nrows, ncols, k+1)
                tsdss.dimshow(img, vmin=-0.1, vmax=1., ticks=False)
            print "saving r-band figure"
            plt.savefig(plotfn)

    return outfns


def stamp_filename(outdir, band, ra, dec):
    return os.path.join(outdir, 'stamp-%s-%.4f-%.4f.fits' % (band, ra, dec))


def make_resampled_psf_images(RCF, band, ra, dec, sdss, 
                            targetwcs, W, H, addToHeader, plots=False,
                            max_exposures=1):
    """ Given a list of (Run, Camcol, Field) tuples, returns a list of
    (img, imgvar, and header) info for stamp sized imgs centered at ra, dec
    """
    # populate list of resampled images and their new psf's
    output_imgs = []

    # zip through each frame, cut out the relevatn patch
    for ifield, (run,camcol,field) in enumerate(RCF[:max_exposures]):
        print """=============================
               RCF %d of %d 
               ======================== """%(ifield, len(RCF))
        # get photofield filename from SDSS, cut it down to relevent RCF
        fn = sdss.retrieve('photoField', run, camcol, field)
        F  = aufits.fits_table(fn)
        F.cut((F.run == run) * (F.camcol == camcol) * (F.field == field))
        print len(F), 'fields'
        assert(len(F) == 1)
        F = F[0]

        # actually get the tractor image (check if it's in cache!)
        boundpixradius = int(np.ceil(np.sqrt(2.) * pixradius))
        print 'RA,Dec,size', (ra, dec, boundpixradius)
        tim, tinfo = tsdss.get_tractor_image_dr9(
            run, camcol, field, band, sdss=sdss, nanomaggies=True,
            roiradecsize=(ra, dec, boundpixradius))
        print 'Got tim:', tim
        frame = sdss.readFrame(run, camcol, field, band)
        if tim is None:
            continue

        # find pixel position for input RA, DEC in tractor image (original field)
        x,y = tim.getWcs().positionToPixel(tsdss.RaDecPos(ra, dec))
        x,y = int(x), int(y)

        # Grab calibration information for header
        tim.sdss_calib   = np.median(frame.getCalibVec())
        tim.sdss_sky     = frame.getSkyAt(x,y)
        iband            = tsdss.band_index(band)
        tim.sdss_gain    = F.gain[iband]
        tim.sdss_darkvar = F.dark_variance[iband]

        # get region of interest in the original frame
        roi = tinfo['roi']
        x0,x1,y0,y1 = roi

        # Resample to common grid
        th,tw = tim.shape
        wwcs = tsdss.TractorWCSWrapper(tim.getWcs(), tw, th)
        try:
            Yo,Xo,Yi,Xi,[rim] = aresample.resample_with_wcs(
                targetwcs, wwcs, [tim.getImage()], Lanczos)
        except aresample.OverlapError:
            continue
        img = np.zeros((H,W))
        img[Yo,Xo] = rim
        iv  = np.zeros((H,W))
        iv[Yo,Xo] = tim.getInvvar()[Yi,Xi]

        # Convert old PSF to new stamp-specific PSF
        newpsf = convert_psf_between_imgs(tim, targetwcs)

        # create the image's header
        hdr = construct_new_header(tim, tinfo, targetwcs, newpsf, 
                                    run, camcol, field, band, addToHeader)

        # add to the list of resampled imgs,
        output_imgs.append((img, iv, hdr))

    return output_imgs


def construct_new_header(tim, tinfo, targetwcs, newpsf, 
                         run, camcol, field, band, addToHeader):
    """ constructs a new header from the old image information (tim, tinfo), 
    writes the new psf, outputs fitsio.FITSHDR object
    """
    hdr = fitsio.FITSHDR()
    targetwcs.add_to_header(hdr)
    hdr.add_record(dict(name='RUN', value=run, comment='SDSS run'))
    hdr.add_record(dict(name='CAMCOL', value=camcol, comment='SDSS camcol'))
    hdr.add_record(dict(name='FIELD', value=field, comment='SDSS field'))
    hdr.add_record(dict(name='BAND', value=band, comment='SDSS band'))

    # Copy from input "frame" header
    orighdr = tinfo['hdr']
    for key in ['NMGY']:
        hdr.add_record(dict(name=key, value=orighdr[key],
                            comment=orighdr.get_comment(key)))
    hdr.add_record(dict(name='CALIB', value=tim.sdss_calib,
                        comment='Mean "calibvec" value for this image'))
    hdr.add_record(dict(name='SKY', value=tim.sdss_sky,
                        comment='SDSS sky estimate at image center'))
    hdr.add_record(dict(name='GAIN', value=tim.sdss_gain,
                        comment='SDSS gain'))
    hdr.add_record(dict(name='DARKVAR', value=tim.sdss_darkvar,
                        comment='SDSS dark variance'))

    # add custom stuff to header
    for (key, value, comment) in addToHeader:
        hdr.add_record(dict(name=key, value=value, comment=comment))

    newpsf.toFitsHeader(hdr, 'PSF_')
    return hdr


def convert_psf_between_imgs(tim, targetwcs):
    """ takes the point spread function MoG from input tractor img (tim), 
    and converts it to a point spread function for the image specified by 
    targetwcs (centered by targetwcs)
    """
    ra, dec  = targetwcs.radec_center()
    th, tw   = tim.shape
    cd       = tim.getWcs().cdAtPixel(tw/2, th/2)
    targetcd = np.array(targetwcs.cd).copy().reshape((2,2))
    trans    = np.dot(np.linalg.inv(targetcd), cd)
    #print 'Tim CD matrix', cd
    #print 'Target CD matrix:', targetcd
    #print 'Transformation matrix:', trans
    psf      = tim.getPsf()
    K        = psf.mog.K
    newmean  = np.zeros_like(psf.mog.mean)
    newvar   = np.zeros_like(psf.mog.var)
    for i,(dx,dy) in enumerate(psf.mog.mean):
        x,y = tim.getWcs().positionToPixel(tsdss.RaDecPos(ra, dec))
        r,d = tim.getWcs().pixelToPosition(x + dx, y + dy)
        #print 'dx,dy', dx,dy
        #print 'ra,dec', r,d
        ok,x0,y0 = targetwcs.radec2pixelxy(ra, dec)
        ok,x1,y1 = targetwcs.radec2pixelxy(r, d)
        #print 'dx2,dy2', x1-x0, y1-y0
        vv = np.array([dx,dy])
        tv = np.dot(trans, vv)
        #print 'dot', tv
        newmean[i,:] = tv

    for i,var in enumerate(psf.mog.var):
        #print 'var', var
        newvar[i,:,:] = np.dot(trans, np.dot(var, trans.T))
        #print 'newvar', newvar[i,:,:]
    newpsf = tsdss.GaussianMixturePSF(psf.mog.amp, newmean, newvar)
    return newpsf


def get_overlapping_run_camcol_field_rerun_301(ra, dec, sdss):
    """ uses window_flist.fits file to find the overlapping frames 
    with the given ra, dec
    """
    wlistfn = sdss.filenames.get('fits_files/window_flist-cut', 'fits_files/window_flist-cut.fits')
    #wfn = os.path.join(os.environ['PHOTO_RESOLVE'], 'window_flist.fits')
    RCF = tsdss.radec_to_sdss_rcf(ra, dec, tablefn=wlistfn)
    print 'Found', len(RCF), 'fields in range.'

    # subselect fields that aren't 157 for some reason...
    keepRCF = []
    for run,camcol,field,r,d in RCF:
        rr = sdss.get_rerun(run, field)
        #print 'Rerun:', rr
        if rr == '157':
            continue
        keepRCF.append((run,camcol,field))
    RCF = keepRCF
    if len(RCF) == 0:
        print 'No run/camcol/fields in rerun 301'
        return
    return RCF


def create_source_table_from_fields(RCF, ra, dec, cutToPrimary,
                                    srcband, sdss):
    # gather rows from sources within a specified radius
    TT = []
    for ifield,(run,camcol,field) in enumerate(RCF):

        # Retrieve SDSS catalog sources in the field
        srcs,objs = tsdss.get_tractor_sources_dr9(
                        run, camcol, field,
                        bandname      = srcband,
                        sdss          = sdss,          # cache is in scratch/
                        radecrad      = (ra, dec, radius*np.sqrt(2.)),
                        nanomaggies   = True,
                        cutToPrimary  = cutToPrimary,
                        getsourceobjs = True,
                        useObjcType   = True)

        print 'Got sources:'
        for src in srcs:
            print '  ', src

        # Write out the sources
        T     = aufits.fits_table()
        T.ra  = [src.getPosition().ra  for src in srcs]
        T.dec = [src.getPosition().dec for src in srcs]

        # same objects, same order
        assert(len(objs) == len(srcs))
        assert(np.all(T.ra == objs.ra))

        # r-band
        bandnum     = 2
        T.primary   = ((objs.resolve_status & 256) > 0)
        T.run       = objs.run
        T.camcol    = objs.camcol
        T.field     = objs.field
        T.is_star   = (objs.objc_type == 6)
        T.frac_dev  = objs.fracdev[:,bandnum]
        T.theta_dev = objs.theta_dev[:,bandnum]
        T.theta_exp = objs.theta_exp[:,bandnum]
        T.phi_dev   = objs.phi_dev_deg[:,bandnum]
        T.phi_exp   = objs.phi_exp_deg[:,bandnum]
        T.ab_dev    = objs.ab_dev[:,bandnum]
        T.ab_exp    = objs.ab_exp[:,bandnum]

        for band in bands:
            bi = tsdss.band_index(band)
            T.set('psfflux_%s' % band, objs.psfflux[:,bi])
            T.set('devflux_%s' % band, objs.devflux[:,bi])
            T.set('expflux_%s' % band, objs.expflux[:,bi])
            T.set('cmodelflux_%s' % band, objs.cmodelflux[:,bi])

        TT.append(T)
    T = tsdss.merge_tables(TT)
    return T


def load_and_clean_objects(object_file, outdir, idx_keep):
    """ loads source information in the object file provided
        inputs:
          - object_file : fits file with columns from PhotoPrimary
          - idx_keep    : index of object to keep

        criteria:
         - removes flagged sources (based on a handful of flags below)
         - psfmag_r < 22.
         - 
    """
    T = aufits.fits_table(object_file)
    print 'Read', len(T), 'objects'
    T.cut(T.nchild == 0)
    print len(T), 'children'
    #T.cut(T.timask == 0)
    #print len(T), 'not in mask'

    #T.cut(np.hypot(T.ra - 5.0562, T.dec - 0.0643) < 0.001)
    # http://skyserver.sdss.org/dr12/en/help/browser/browser.aspx#&&history=enum+PhotoFlags+E
    for flagname,flagval in [('BRIGHT', 0x2),
                             ('EDGE', 0x4),
                             ('NODEBLEND', 0x40),
                             ('DEBLEND_TOO_MANY_PEAKS' , 0x800),
                             ('NOTCHECKED', 0x80000),
                             ('TOO_LARGE', 0x1000000),
                             ('BINNED2', 0x20000000),
                             ('BINNED4', 0x40000000),
                             ('SATUR_CENTER', 0x80000000000),
                             ('INTERP_CENTER', 0x100000000000),
                             ('MAYBE_CR', 0x100000000000000),
                             ('MAYBE_EGHOST', 0x200000000000000),
                  ]:
        T.cut(T.flags & flagval == 0)
        print len(T), 'without', flagname, 'bit set'
        #pass

    # Cut to objects that are likely to appear in the individual images
    #T.cut(T.psfmag_r < 22.)
    #print 'Cut to', len(T), 'with psfmag_r < 22 in coadd'

    # construct labels, write to a single stamps.fits file
    idx_keep = int(idx_keep)
    assert idx_keep <= len(T), "only %d sources in final cut, idx_keep is too big"%len(T)
    T.tag = np.array(['%.4f-%.4f.fits' % (r,d) for r,d in zip(T.ra, T.dec)])
    #T[:Nkeep].writeto(os.path.join(outdir, 'stamps.fits'),
    #                  columns='''tag objid run camcol field ra dec
    #                             psfmag_u psfmag_g psfmag_r psfmag_i psfmag_z
    #                             modelmag_u modelmag_g modelmag_r
    #                             modelmag_i modelmag_z'''.split())
    return T[idx_keep:idx_keep+1]

def write_cat_files(T, outdir):
    """ writes out single source files for each object in T, to outdir
    """
    # Write out Stripe82 measurements...
    radius = np.sqrt(2.) * pixradius * pixscale / 3600.

    for i in range(len(T)):
        # looks for sources nearby T[i], within radius.  it returns index in 
        # both first (short) list and second (long) list
        I,J,d = asphere.match_radec(np.array([T.ra[i]]), np.array([T.dec[i]]),
                                              T.ra, T.dec, radius)
        print len(J), 'matched within', radius*3600., 'arcsec'
        t = T[J]
        print len(t), 'matched within', radius*3600., 'arcsec'

        tt = aufits.fits_table()
        cols = ['ra','dec','run','camcol','field',#'probpsf',
                #'flags', #'type',
                'fracdev_r', #'probpsf_r', 
                'devrad_r','devraderr_r', 'devab_r', 'devaberr_r',
                'devphi_r', 'devphierr_r',
                'exprad_r','expraderr_r', 'expab_r', 'expaberr_r',
                'expphi_r', 'expphierr_r',
                ]
        for c in cols:
            cout = c
            # drop "_r" from dev/exp shapes
            if cout.endswith('_r'):
                cout = cout[:-2]

            coutmap = dict(devrad='theta_dev',
                           devphi='phi_dev',
                           devab ='ab_dev',
                           devraderr='theta_dev_err',
                           devphierr='phi_dev_err',
                           devaberr ='ab_dev_err',
                           exprad='theta_exp',
                           expphi='phi_exp',
                           expab ='ab_exp',
                           expraderr='theta_exp_err',
                           expphierr='phi_exp_err',
                           expaberr ='ab_exp_err',
                           fracdev='frac_dev')
            cout = coutmap.get(cout, cout)
            tt.set(cout, t.get(c))

        tt.is_star = (t.type == 6)

        for magname in ['psf', 'dev', 'exp']:
            for band in 'ugriz':
                mag    = t.get('%smag_%s' % (magname, band))
                magerr = t.get('%smagerr_%s' % (magname, band))

                ### FIXME -- arcsinh mags??
                flux  = tsdss.NanoMaggies.magToNanomaggies(mag)
                dflux = np.abs(flux * np.log(10.)/-2.5 * magerr)

                tt.set('%sflux_%s' % (magname, band), flux)
                tt.set('%sfluxerr_%s' % (magname, band), dflux)

        for band in 'ugriz':
            # http://www.sdss3.org/dr10/algorithms/magnitudes.php#cmodel
            fexp = tt.get('expflux_%s' % band)
            fdev = tt.get('expflux_%s' % band)
            fracdev = t.get('fracdev_%s' % band)
            tt.set('cmodelflux_%s' % band, fracdev * fdev + (1.-fracdev) * fexp)

        catfn = os.path.join(outdir, 'cat-s82-%.4f-%.4f.fits' % (t.ra[0], t.dec[0]))
        tt.writeto(catfn)
        print 'Wrote', catfn


if __name__ == '__main__':
    main()



#if plots:
#    plt.clf()
#    img = tim.getImage()
#    mn,mx = [np.percentile(img,p) for p in [25,99]]
#    tsdss.dimshow(img, vmin=mn, vmax=mx)
#    xx,yy = [],[]
#    for src in srcs:
#        x,y = tim.getWcs().positionToPixel(src.getPosition())
#        xx.append(x)
#        yy.append(y)
#    ax = plt.axis()
#    plt.plot(xx, yy, 'r+')
#    plt.axis(ax)
#    plt.savefig('tim-%s%i.png' % (band, ifield))


#if plots:
#    plt.clf()
#    mn,mx = [np.percentile(img,p) for p in [25,99]]
#    tsdss.dimshow(img, vmin=mn, vmax=mx)
#    xx,yy = [],[]
#    for src in srcs:
#        rd = src.getPosition()
#        ok,x,y = targetwcs.radec2pixelxy(rd.ra, rd.dec)
#        xx.append(x-1)
#        yy.append(y-1)
#    ax = plt.axis()
#    plt.plot(xx, yy, 'r+')
#    plt.axis(ax)
#    plt.savefig('rim-%s%i.png' % (band, ifield))


