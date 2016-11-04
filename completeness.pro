pro completeness

; End-to-end example of a CPROPSTOO run on an unnamed ALMA data set.

; -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
; PREP
; -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

; The cube needs to have km/s velocity units with units of Kelvin in
; order to proceed. "prep_cube" is a necessarily imperfect attempt to
; clean up a header to CPROP standrds. Often the user will need to do
; a bit more work at this stage. Then "cprops_check_header" will
; attempt to evaluate the validity of the header for further
; analysis. Note that this stage is not strictly necessary and,
; because of the wide variety of possible headers, may often fail for
; headers far off what is expected.

;  rootlist = ['m83.fakes_5.00', $
;              'm83.fakes_6.00'] ;, $
;              'm83.fakes_7.0']
  rootlist = 'm83.fakes_'+string(indgen(51)/50.+5.0,format='(F4.2)')

  ;; rootlist = ['m83.fakes_4.0',$
  ;;             'm83.fakes_4.1',$
  ;;             'm83.fakes_4.2',$
  ;;             'm83.fakes_4.3',$
  ;;             'm83.fakes_4.4',$
  ;;             'm83.fakes_4.5',$
  ;;             'm83.fakes_4.6',$
  ;;             'm83.fakes_4.7',$
  ;;             'm83.fakes_4.8',$
  ;;             'm83.fakes_4.9',$
  ;;             'm83.fakes_5.0',$
  ;;             'm83.fakes_5.1',$
  ;;             'm83.fakes_5.2',$
  ;;             'm83.fakes_5.3',$
  ;;             'm83.fakes_5.4',$
  ;;             'm83.fakes_5.5',$
  ;;             'm83.fakes_5.6',$
  ;;             'm83.fakes_5.7',$
  ;;             'm83.fakes_5.8',$
  ;;             'm83.fakes_5.9',$
  ;;             'm83.fakes_6.0',$
  ;;             'm83.fakes_6.1',$
  ;;             'm83.fakes_6.2',$
  ;;             'm83.fakes_6.3',$
  ;;             'm83.fakes_6.4',$
  ;;             'm83.fakes_6.5']
  ;;             'm83.fakes_6.6',$
  ;;             'm83.fakes_6.7',$
  ;;             'm83.fakes_6.8',$
  ;;             'm83.fakes_6.9',$
  ;;             'm83.fakes_7.0']
            
  rootlist = reverse(rootlist)
  for i = 0,n_elements(rootlist)-1 do begin
     
     root = rootlist[i]

; CHECK THE HEADER
     print, cprops_check_header(in_file="../data/"+root+".fits" $
                                , perfect=perfect $
                                , comments=comments)

; PRINT ANY COMMENTS FROM THE HEADER CHECK
     print, comments

; MAKE A PRIMARY-BEAM CORRECTED CUBE. WE WILL USE THIS FOR
; MEASUREMENTS BUT NOT FOR MASKING OR ASSIGNMENT (FLATTER NOISE MAKES
; THOSE EASIER).

     cube = readfits("../data/"+root+".fits", hdr)
;  flux = readfits("../data/"+root+"_flux.fits", hdr)
;  writefits, "../data/"+root+"_pbcorr.fits", cube/flux, hdr

; -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
; MAKE A MASK
; -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

; Estimate the noise from the data. We will use this to make a mask
; identifying bright emission. The noise estimate can have a wide
; range of complexities. Here the call uses a moving box of 9x9 pixels
; and calculates a two-d noise map. Optionally, it can fit a spectral
; dependence of the noise or just treat the whole cube as a single
; distribution. The example case here would probably be better served
; with a single value.

     ;; make_noise_cube $
     ;;    , cube_file = "../data/"+root+".fits" $
     ;;    , out_file = "../data/"+root+"_noise.fits" $
     ;;    , box = 9 $
     ;;    , /iterate $
     ;;    , /twod_only $
     ;;    , /collapse

; Now make a mask of bright emission. Require two channels above
; signal-to-noise five and expand these regions into a surface with at
; least two channels above signal-to-noise two.

     make_cprops_mask $
        , infile = "../data/"+root+".fits" $
        , outfile="../data/"+root+"_mask.fits" $   
        , rmsfile = "../data/m83.fakes.noise.fits" $ ;"+root+"_noise.fits" $
        , hi_thresh = 2.5 $
        , hi_nchan = 2 $
        , lo_thresh = 1.0 $
        , lo_nchan = 1 $
        , min_pix = 40

; Rerun the masking with the cube "inverted" (flipped to negative) in
; order to check the stringency of the mask. The number of pixels in
; the mask here give an estimate of the flase positive rate.

     make_cprops_mask $
        , infile = "../data/"+root+".fits" $
        , outfile="../data/"+root+"_falsepos_test.fits" $   
        , rmsfile = "../data/m83.fakes.noise.fits" $ ;"+root+"_noise.fits" $
        , hi_thresh = 2.5 $
        , hi_nchan = 2 $
        , lo_thresh = 1.0 $
        , lo_nchan = 1 $
        , /invert $
        , min_pix = 40

; ... for our case we can load this cube and see that the mask is
; empty. To some approximation we expect no false positives.  

; -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
; IDENTIFY LOCAL MAXIMA IN THE CUBE
; -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

; Identify local maxima in the cube. The program slides a sub-cube
; through the masked region and identifies significant local
; maxima. The results are saved in either a text file or IDL save file
; for further future analysis.
;
; This step allows a good deal of user input in what defines a
; significant local maximum. DELTA defines the contrast required above
; any contour shared with another maximum (this can be defined in
; intensity or SNR units). MINPIX and MINAREA define the area and
; volume required to be uniquely associated with an individual maximum
; for it to be valid. SPECFRIENDS and FRIENDS define the search box in
; pixel units. Think about this as having knobs to turn the contrast
; with the background, the search area, and the region required for a
; maximum to be valid.

     find_local_max $
        , infile="../data/"+root+".fits" $
        , inmask ="../data/"+root+"_mask.fits" $
        , text_out="../measurements/"+root+"_lmax.txt" $
        , /verbose $
        , delta=2.5 $
        , /snr $
        , minpix=40 $
        , minarea=8 $
        , minvchan=2 $
        , friends=11 $
        , specfriends=2

; Visualize the selected peaks:

     ;; peak_map = max(readfits("../data/"+root+".fits", hdr),dim=3,/nan)
     ;; loadct, 33
     ;; disp, peak_map, /sq
     ;; readcol, "../measurements/"+root+"_lmax.txt", comment="#", x, y
     ;; loadct, 0
     ;; oplot, x, y, color=255, ps=1

; Manipulate the FIND_LOCAL_MAX call and iterate until happy.

; -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
; ASSIGN VOXELS IN THE CUBE TO LOCAL MAXIMA
; -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
     
; Based on the list of local maxima and the mask, "color" the cube
; with assignments of voxels to discrete objects. The result is a cube
; with pixels corresponding to object number. The current propgram
; includes two methodologies to do this: the "CPROPS" way of assigning
; only emission uniquely associated with a maximum (defined as
; emission above any shared contour with another maximum) and the
; "CLUMPFIND" way of assigning emission within a shared contour to the
; nearest maximum on the basis of proximity.

;; ; CPROPS ASSIGNMENT
;;   assign_cprops $
;;      , kernfile="../measurements/"+root+"_lmax.txt" $
;;      , infile="../data/"+root+".fits" $
;;      , inmask="../data/"+root+"_mask.fits" $
;;      , outfile="../data/"+root+"_assign_cprops.fits" $
;;      , /verbose

; CLUMPFIND ASSIGNMENT
     assign_clfind $
        , kernfile="../measurements/"+root+"_lmax.txt" $
        , infile="../data/"+root+".fits" $
        , inmask="../data/"+root+"_mask.fits" $
        , outfile="../data/"+root+"_assign_clfind.fits" $
        , /verbose

; -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
; MEASURE CLOUD MOMENTS
; -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
     
; Using the data cube (note that we now use the primary-beam corrected
; cube) and the assignment cubes generated in the last step measure
; moments for each object in the assignment cube. Here moments is used
; loosely to refer to pixelwise properties (rather than properties
; with physical units, which we compute in the next step). 
;
; Additional analysis code could also be plugged in here to operate on
; the assignment+data combination.

;; ; THE CPROPS CASE
;;   cube_to_moments $
;;      , infile="../data/"+root+".fits" $
;;      , inassign="../data/"+root+"_assign_cprops.fits" $
;;      , outfile="../measurements/"+root+"_moments_cprops.idl" $
;;      , rmsfile = "../data/"+root+"_noise.fits" $
;;      , /verbose           

; THE CLUMPFIND CASE
     cube_to_moments $
        , infile="../data/"+root+".fits" $
        , inassign="../data/"+root+"_assign_clfind.fits" $
        , outfile="../measurements/"+root+"_moments_clfind.idl" $
        , rmsfile = "../data/m83.fakes.noise.fits" $ ;"+root+"_noise.fits" $
        , /verbose           

; -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
; CONVERT MOMENTS TO PROPERTIES
; -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

; We now have a structure of moments for each object. The last step
; before science is to convert this into physical units. This step
; accepts moment measurements and uses a header and other assumptions
; (e.g., distance, conversion factor, etc.) to write a "properties"
; IDL structure. We run it once each for the two assignment cubes.

;; ; THE CPROPS CASE
;;   moments_to_props $
;;      , hdrfile="../data/"+root+".fits" $
;;      , infile="../measurements/"+root+"_moments_cprops.idl" $
;;      , fits_file="../measurements/"+root+"_props_cprops.fits" $
;;      , dist=4.6d6 $
;;      , /verbose

; THE CLUMPFIND CASE
     moments_to_props $
        , hdrfile="../data/"+root+".fits" $
        , infile="../measurements/"+root+"_moments_clfind.idl" $
        , fits_file="../measurements/"+root+"_props_clfind.fits" $
        , dist=4.6d6 $
        , /verbose

;; ; -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
;; ; EXTRACT PROPERTIES FOR ANOTHER LINE
;; ; -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

;; ; A major advantage of this modular approach is the ability to apply
;; ; assignments to a data cube for another line.  Here we provide an
;; ; example call for a 13CO (3-2) data cube. Note that we require
;; ; matched astrometry in all three axes for the assignment cube to be
;; ; valid. Alignment options exist within CASA, IDL, etc.
     
;; ; Notice that the cube_to_moments call uses the data file for the 13CO
;; ; line but the assignment from the 12CO line.

;; ; PRIMARY BEAM CORRECTED 13CO CUBE
;;   cube = readfits("../data/13"+root+".fits", hdr)
;;   flux = readfits("../data/13"+root+"_flux.fits", hdr)
;;   writefits, "../data/13"+root+".fits", cube/flux, hdr

;; ; CALCULATE CPROPS ASSIGNMENT MOMENTS FOR 13CO
;;   cube_to_moments $
;;      , infile="../data/13"+root+".fits" $
;;      , inassign="../data/"+root+"_assign_cprops.fits" $
;;      , outfile="../measurements/13"+root+"_moments_cprops.idl" $
;;      , /verbose           

;; ; CALCULATE CLUMPFIND ASSIGNMENT MOMENTS FOR 13CO
;;   cube_to_moments $
;;      , infile="../data/13"+root+".fits" $
;;      , inassign="../data/"+root+"_assign_clfind.fits" $
;;      , outfile="../measurements/13"+root+"_moments_clfind.idl" $
;;      , /verbose           

;; ; CALCULATE CPROPS PROPERTIES FOR 13CO
;;   moments_to_props $
;;      , hdrfile="../data/13"+root+".fits" $
;;      , infile="../measurements/13"+root+"_moments_cprops.idl" $
;;      , idl_file="../measurements/13"+root+"_props_cprops.idl" $
;;      , dist=11.9d6 $
;;      , /verbose

;; ; CALCULATE CLUMPFIND PROPERTIES FOR 13CO
;;   moments_to_props $
;;      , hdrfile="../data/13"+root+".fits" $
;;      , infile="../measurements/13"+root+"_moments_clfind.idl" $
;;      , idl_file="../measurements/13"+root+"_props_clfind.idl" $
;;      , dist=11.9d6 $
;;      , /verbose

; -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
; DENDROGRAMS
; -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

; The code also contains the machinery for an approach based on the
; "dendrograms" formalism. In this case moments are measured for each
; combination of local maximum and contour level and then saved. These
; can be analyzed or visualized using the dendrogram tree diagrams or
; some other approach.

;; ; GET THE MOMENTS BY LEVEL
;;      cube_to_level_moments $
;;         , kernfile = "../measurements/"+root+"_lmax.txt" $
;;         , infile="../data/"+root+".fits" $
;;         , inmask="../data/"+root+"_mask.fits" $
;;         , outfile="../measurements/"+root+"_level_moments.idl" $
;;         , /verbose
     
;; ; GET THE PROPERTIES FOR THOSE MOMENTS  
;;      moments_to_props $
;;         , hdrfile="../data/"+root+".fits" $
;;         , infile="../measurements/"+root+"_level_moments.idl" $
;;         , idl_file="../measurements/"+root+"_level_props.idl" $
;;         , dist=4.6d6 $
;;         , /verbose

; -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
; ... PROFIT?
; -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

; You now have IDL save files holding the property measurements for
; the objects extraced from your data set. You can take these and use
; them to do science.
  endfor
end
