    !     Code for Anisotropies in the Microwave Background
    !     by Antony Lewis (http://cosmologist.info/) and Anthony Challinor
    !     See readme.html for documentation. This is a sample driver routine that reads
    !     in one set of parameters and produdes the corresponding output.

    program driver
    use IniFile
    use CAMB
    use LambdaGeneral
    use Lensing
    use AMLUtils
    use Transfer
    use constants
    use Bispectrum
    use CAMBmain
    use NonLinear
    
!#SimDataADD
    use CSGalCalc
    use Clusters
!#SimDataAdd
    
#ifdef NAGF95
    use F90_UNIX
#endif
    implicit none

    Type(CAMBparams) P

    character(LEN=Ini_max_string_len) numstr, VectorFileName, &
        InputFile, ScalarFileName, TensorFileName, TotalFileName, LensedFileName,&
        LensedTotFileName, LensPotentialFileName,ScalarCovFileName
    integer i
    character(LEN=Ini_max_string_len) TransferFileNames(max_transfer_redshifts), &
        MatterPowerFileNames(max_transfer_redshifts), outroot, version_check
    real(dl) output_factor, nmassive

#ifdef WRITE_FITS
    character(LEN=Ini_max_string_len) FITSfilename
#endif

    logical bad
!#SimDataAdd
    logical LSS_CL_Differ
!#SimDataAdd

    InputFile = ''
    if (GetParamCount() /= 0)  InputFile = GetParam(1)
    if (InputFile == '') error stop 'No parameter input file'

    call Ini_Open(InputFile, 1, bad, .false.)
    if (bad) error stop 'Error opening parameter file'

    Ini_fail_on_not_found = .false.

    outroot = Ini_Read_String('output_root')
    if (outroot /= '') outroot = trim(outroot) // '_'

    highL_unlensed_cl_template = Ini_Read_String_Default('highL_unlensed_cl_template',highL_unlensed_cl_template)

    call CAMB_SetDefParams(P)

    P%WantScalars = Ini_Read_Logical('get_scalar_cls')
    P%WantVectors = Ini_Read_Logical('get_vector_cls',.false.)
    P%WantTensors = Ini_Read_Logical('get_tensor_cls',.false.)

    P%OutputNormalization=outNone
    output_factor = Ini_Read_Double('CMB_outputscale',1.d0)

    P%WantCls= P%WantScalars .or. P%WantTensors .or. P%WantVectors

    P%PK_WantTransfer=Ini_Read_Logical('get_transfer')

    AccuracyBoost  = Ini_Read_Double('accuracy_boost',AccuracyBoost)
    lAccuracyBoost = Ini_Read_Real('l_accuracy_boost',lAccuracyBoost)
    HighAccuracyDefault = Ini_Read_Logical('high_accuracy_default',HighAccuracyDefault)

    P%NonLinear = Ini_Read_Int('do_nonlinear',NonLinear_none)

    P%DoLensing = .false.
    if (P%WantCls) then
        if (P%WantScalars  .or. P%WantVectors) then
            P%Max_l = Ini_Read_Int('l_max_scalar')
            P%Max_eta_k = Ini_Read_Double('k_eta_max_scalar',P%Max_l*2._dl)
            if (P%WantScalars) then
                P%DoLensing = Ini_Read_Logical('do_lensing',.false.)
                if (P%DoLensing) lensing_method = Ini_Read_Int('lensing_method',1)
            end if
            if (P%WantVectors) then
                if (P%WantScalars .or. P%WantTensors) error stop 'Must generate vector modes on their own'
                i = Ini_Read_Int('vector_mode')
                if (i==0) then
                    vec_sig0 = 1
                    Magnetic = 0
                else if (i==1) then
                    Magnetic = -1
                    vec_sig0 = 0
                else
                    error stop 'vector_mode must be 0 (regular) or 1 (magnetic)'
                end if
            end if
        end if

        if (P%WantTensors) then
            P%Max_l_tensor = Ini_Read_Int('l_max_tensor')
            P%Max_eta_k_tensor =  Ini_Read_Double('k_eta_max_tensor',Max(500._dl,P%Max_l_tensor*2._dl))
        end if
    endif
    
   !#SimDataADD
    P%DoCS = Ini_Read_Logical('do_cosmic_shear',.false.)
    P%DoGal = Ini_Read_Logical('do_galaxy_ps',.false.)
    P%DoCSXGal = Ini_Read_Logical('do_shear_x_galaxy_ps',.false.)
    P%DoClusters = Ini_Read_Logical('do_clusters',.false.) 
    P%DoCMB = Ini_Read_Logical('do_cmb',.false.)
    LSS_CL_Differ = Ini_Read_Logical('lss_and_cl_differ',.false.)
    P%OutputSimDataFiles = Ini_Read_Logical('output_sim_data_files',.false.)
    P%sim_random_cmb = Ini_Read_Logical('sim_random_cmb',.false.)
    P%sim_random_seed = Ini_Read_Int('sim_random_seed')
    setname = Ini_Read_String('output_root')
    datarep = Ini_Read_String('datarep')
    

    if (P%DoCS .or. P%DoGal .or. P%DoClusters) P%DoShePowFoC = .true.

    if (P%DoCSXGal .and. .not. (P%DoCS .and. P%DoGal)) then
        stop 'Sorry, cross correlation requires shear and galaxy power spectrum calculations'
    end if

    if(P%DoShePowFoC .and. .not. P%WantScalars) then
        stop 'Sorry, ShePowFoC requires WantScalars = T'
    end if

    if(P%DoShePowFoC) then
        if(LSS_CL_Differ) then
            P%ngal_lss = Ini_Read_Double('ngal_lss')
            P%nz_alpha_lss = Ini_Read_Double('nz_alpha_lss')
            P%nz_beta_lss = Ini_Read_Double('nz_beta_lss')
            P%nz_z0_lss = Ini_Read_Double('nz_z0_lss')
            P%mean_int_ellip_lss = Ini_Read_Double('mean_int_ellip_lss')
            P%fskylss = Ini_Read_Double('fsky_lss')
            P%ngal_cl = Ini_Read_Double('ngal_cl')
            P%nz_alpha_cl = Ini_Read_Double('nz_alpha_cl')
            P%nz_beta_cl = Ini_Read_Double('nz_beta_cl')
            P%nz_z0_cl = Ini_Read_Double('nz_z0_cl')
            P%mean_int_ellip_cl = Ini_Read_Double('mean_int_ellip_cl')
            P%fskycl = Ini_Read_Double('fsky_cl')
        else
            P%ngal_lss = Ini_Read_Double('ngal')
            P%nz_alpha_lss = Ini_Read_Double('nz_alpha')
            P%nz_beta_lss = Ini_Read_Double('nz_beta')
            P%nz_z0_lss = Ini_Read_Double('nz_z0')
            P%mean_int_ellip_lss = Ini_Read_Double('mean_int_ellip')
            P%fskylss = Ini_Read_Double('fsky')
            P%ngal_cl = Ini_Read_Double('ngal')
            P%nz_alpha_cl = Ini_Read_Double('nz_alpha')
            P%nz_beta_cl = Ini_Read_Double('nz_beta')
            P%nz_z0_cl = Ini_Read_Double('nz_z0')
            P%mean_int_ellip_cl = Ini_Read_Double('mean_int_ellip')
            P%fskycl = Ini_Read_Double('fsky')
        end if
        P%photo_error = Ini_Read_Double('photo_error',0.d0)
    end if

    if(P%DoCMB) then
        P%lmin_cmb = Ini_Read_Int('lmin_CMB')
        P%lmax_cmb = Ini_Read_Int('lmax_CMB')
        P%fskycmb = Ini_Read_Double('fsky_cmb')
        P%sim_ncmbcls = Ini_Read_Int('ncmbcls')
	P%nchan = Ini_Read_Int('nchan')

        allocate(P%fwhm_arcmin(1:P%nchan),P%Noise_Sigma_T(1:P%nchan),P%Noise_Sigma_P(1:P%nchan))

        numstr = Ini_Read_String('fwhm_arcmin')
        read(numstr,*) P%fwhm_arcmin
        numstr = Ini_Read_String('sigma_T')
        read(numstr,*) P%Noise_Sigma_T
        numstr = Ini_Read_String('sigma_P')
        read(numstr,*) P%Noise_Sigma_P
    end if

    if(P%DoCS) then
        P%nTomoBin(1) = Ini_Read_Int('n_tomo_bin_shear',0)
        if (P%nTomoBin(1) .eq. 0) stop 'Please specify the number of shear tomography bins'

        auto_binning(1) = Ini_Read_Logical('auto_binning_shear')
        if (auto_binning(1)) then
             auto_binning_zmin(1) = Ini_Read_Double('auto_binning_zmin_shear')
             auto_binning_zmax(1) = Ini_Read_Double('auto_binning_zmax_shear')
        else
             numstr = Ini_Read_String('zph_low_shear')
             read(numstr,*) P%zph_low(1:P%nTomoBin(1),1)
             numstr=Ini_Read_String('zph_high_shear')
             read(numstr,*) P%zph_high(1:P%nTomoBin(1),1)
        end if

        P%lmin_CS = Ini_Read_Int('lmin_CS',2)
        P%lmax_CS = Ini_Read_Int('lmax_CS',0)
        if(P%lmax_CS .eq. 0) stop 'Please specify maximum multipole for shear power spectra'
    end if

    if(P%DoGal) then
        P%nTomoBin(2)=Ini_Read_Int('n_tomo_bin_galPS',0)
        if (P%nTomoBin(2) .eq. 0) stop 'Please specify the number of galaxy power spectrum&
           & tomography bins'

        auto_binning(2) = Ini_Read_Logical('auto_binning_galPS')
        if (auto_binning(2)) then
            auto_binning_zmin(2) = Ini_Read_Double('auto_binning_zmin_galPS')
            auto_binning_zmax(2) = Ini_Read_Double('auto_binning_zmax_galPS')
        else
            numstr=Ini_Read_String('zph_low_galPS')
            read(numstr,*) P%zph_low(1:P%nTomoBin(2),2)
            numstr=Ini_Read_String('zph_high_galPS')
            read(numstr,*) P%zph_high(1:P%nTomoBin(2),2)
        end if

        P%lmin_gal = Ini_Read_Int('lmin_gal',2)
        P%lmax_gal = Ini_Read_Int('lmax_gal',0)
        if(P%lmax_gal .eq. 0) stop 'Please specify maximum multipole for galaxy power spectra'
    end if

    if(P%DoClusters) then 
        P%cl_zbins = Ini_Read_Int('n_tomo_bin_cl')
        P%cl_Mbins = Ini_Read_Int('n_mass_bin_cl')
        P%c_nfw = Ini_Read_Double('c_nfw')
        P%theta_G = Ini_Read_Double('theta_G')
        P%lnM_sigma = Ini_Read_Double('lnM_sigma')
        P%signal2noise = Ini_Read_Double('signal2noise')
        P%completeness = Ini_Read_Double('completeness')
        P%efficiency = Ini_Read_Double('efficiency')
        P%zmin_cl = Ini_Read_Double('zmin_cl')
        P%zmax_cl = Ini_Read_Double('zmax_cl')
    end if
!#SimDataADD
 

    !  Read initial parameters.

!#SimDataReplace
    !call DarkEnergy_ReadParams(DefIni)
    P%w0 = Ini_Read_Double('w0')
    P%wa = Ini_Read_Double('wa')
    P%cs2_de = Ini_Read_Double('cs2_de')
!#SimDataReplace
    
    P%h0     = Ini_Read_Double('hubble')

    if (Ini_Read_Logical('use_physical',.false.)) then
        P%omegab = Ini_Read_Double('ombh2')/(P%H0/100)**2
        P%omegac = Ini_Read_Double('omch2')/(P%H0/100)**2
        P%omegan = Ini_Read_Double('omnuh2')/(P%H0/100)**2
        P%omegav = 1- Ini_Read_Double('omk') - P%omegab-P%omegac - P%omegan
    else
        P%omegab = Ini_Read_Double('omega_baryon')
        P%omegac = Ini_Read_Double('omega_cdm')
        P%omegav = Ini_Read_Double('omega_lambda')
        P%omegan = Ini_Read_Double('omega_neutrino')
    end if

    P%tcmb   = Ini_Read_Double('temp_cmb',COBE_CMBTemp)
    P%yhe    = Ini_Read_Double('helium_fraction',0.24_dl)
    P%Num_Nu_massless  = Ini_Read_Double('massless_neutrinos')

    P%Nu_mass_eigenstates = Ini_Read_Int('nu_mass_eigenstates',1)
    if (P%Nu_mass_eigenstates > max_nu) error stop 'too many mass eigenstates'

    numstr = Ini_Read_String('massive_neutrinos')
    read(numstr, *) nmassive
    if (abs(nmassive-nint(nmassive))>1e-6) error stop 'massive_neutrinos should now be integer (or integer array)'
    read(numstr,*, end=100, err=100) P%Nu_Mass_numbers(1:P%Nu_mass_eigenstates)
    P%Num_Nu_massive = sum(P%Nu_Mass_numbers(1:P%Nu_mass_eigenstates))

    if (P%Num_Nu_massive>0) then
        P%share_delta_neff = Ini_Read_Logical('share_delta_neff', .true.)
        numstr = Ini_Read_String('nu_mass_degeneracies')
        if (P%share_delta_neff) then
            if (numstr/='') write (*,*) 'WARNING: nu_mass_degeneracies ignored when share_delta_neff'
        else
            if (numstr=='') error stop 'must give degeneracies for each eigenstate if share_delta_neff=F'
            read(numstr,*) P%Nu_mass_degeneracies(1:P%Nu_mass_eigenstates)
        end if
        numstr = Ini_Read_String('nu_mass_fractions')
        if (numstr=='') then
            if (P%Nu_mass_eigenstates >1) error stop 'must give nu_mass_fractions for the eigenstates'
            P%Nu_mass_fractions(1)=1
        else
            read(numstr,*) P%Nu_mass_fractions(1:P%Nu_mass_eigenstates)
        end if
    end if
    
!#SimDataAdd
    if(P%DoShePowFoC) then
        call SetTomoBins(P) !Set tomography bins
        write (*,*) 'overriding transfer settings for Cosmic Shear calcs'
        P%WantTransfer = .true.
        call Transfer_SetForNonlinearLensingCS(P%Transfer)
        P%Transfer%PK_redshifts = P%Transfer%redshifts
        P%Transfer%PK_num_redshifts = P%Transfer%num_redshifts
        P%Transfer%high_precision = Ini_Read_Logical('transfer_high_precision',.false.)
        call Transfer_SortAndIndexRedshifts(P%Transfer)
    else
!#SimDataAdd

    !JD 08/13 begin changes for nonlinear lensing of CMB + LSS compatibility
    !P%Transfer%redshifts -> P%Transfer%PK_redshifts and P%Transfer%num_redshifts -> P%Transfer%PK_num_redshifts
    !in the P%WantTransfer loop.
    
    if (((P%NonLinear==NonLinear_lens .or. P%NonLinear==NonLinear_both) .and. P%DoLensing) .or. P%PK_WantTransfer) then
        P%Transfer%high_precision=  Ini_Read_Logical('transfer_high_precision',.false.)
    else
        P%transfer%high_precision = .false.
    endif
    if (P%PK_WantTransfer) then
        P%Transfer%accurate_massive_neutrinos = Ini_Read_Logical('accurate_massive_neutrino_transfers',.false.)
    else
        P%Transfer%accurate_massive_neutrinos = .false.
    end if

    if (P%NonLinear/=NonLinear_none) call NonLinear_ReadParams(DefIni)

    if (P%PK_WantTransfer)  then
        P%WantTransfer  = .true.
        P%transfer%kmax          =  Ini_Read_Double('transfer_kmax')
        P%transfer%k_per_logint  =  Ini_Read_Int('transfer_k_per_logint')
        P%transfer%PK_num_redshifts =  Ini_Read_Int('transfer_num_redshifts')

        transfer_interp_matterpower = Ini_Read_Logical('transfer_interp_matterpower', transfer_interp_matterpower)
        transfer_power_var = Ini_read_int('transfer_power_var',transfer_power_var)
        if (P%transfer%PK_num_redshifts > max_transfer_redshifts) error stop 'Too many redshifts'
        do i=1, P%transfer%PK_num_redshifts
            P%transfer%PK_redshifts(i)  = Ini_Read_Double_Array('transfer_redshift',i,0._dl)
            transferFileNames(i)     = Ini_Read_String_Array('transfer_filename',i)
            MatterPowerFilenames(i)  = Ini_Read_String_Array('transfer_matterpower',i)
            if (TransferFileNames(i) == '') then
                TransferFileNames(i) =  trim(numcat('transfer_',i))//'.dat'
            end if
            if (MatterPowerFilenames(i) == '') then
                MatterPowerFilenames(i) =  trim(numcat('matterpower_',i))//'.dat'
            end if
            if (TransferFileNames(i)/= '') &
                TransferFileNames(i) = trim(outroot)//TransferFileNames(i)
            if (MatterPowerFilenames(i) /= '') &
                MatterPowerFilenames(i)=trim(outroot)//MatterPowerFilenames(i)
        end do
    else
        P%Transfer%PK_num_redshifts = 1
        P%Transfer%PK_redshifts = 0
    end if

    if ((P%NonLinear==NonLinear_lens .or. P%NonLinear==NonLinear_both) .and. P%DoLensing) then
        P%WantTransfer  = .true.
        call Transfer_SetForNonlinearLensing(P%Transfer)
    end if

    call Transfer_SortAndIndexRedshifts(P%Transfer)
    !JD 08/13 end changes

    P%transfer%kmax=P%transfer%kmax*(P%h0/100._dl)
    
!#SimDataAdd
    end if
!#SimDataAdd

    Ini_fail_on_not_found = .false.

    DebugParam = Ini_Read_Double('DebugParam',DebugParam)
    ALens = Ini_Read_Double('Alens',Alens)

    call Reionization_ReadParams(P%Reion, DefIni)
    call InitialPower_ReadParams(P%InitPower, DefIni, P%WantTensors)
    call Recombination_ReadParams(P%Recomb, DefIni)
    if (Ini_HasKey('recombination')) then
        i = Ini_Read_Int('recombination',1)
        if (i/=1) error stop 'recombination option deprecated'
    end if

    call Bispectrum_ReadParams(BispectrumParams, DefIni, outroot)

    if (P%WantScalars .or. P%WantTransfer) then
        P%Scalar_initial_condition = Ini_Read_Int('initial_condition',initial_adiabatic)
        if (P%Scalar_initial_condition == initial_vector) then
            P%InitialConditionVector=0
            numstr = Ini_Read_String('initial_vector',.true.)
            read (numstr,*) P%InitialConditionVector(1:initial_iso_neutrino_vel)
        end if
        if (P%Scalar_initial_condition/= initial_adiabatic) use_spline_template = .false.
    end if

    if (P%WantScalars) then
        ScalarFileName = trim(outroot)//Ini_Read_String('scalar_output_file')
        LensedFileName =  trim(outroot) //Ini_Read_String('lensed_output_file')
        LensPotentialFileName =  Ini_Read_String('lens_potential_output_file')
        if (LensPotentialFileName/='') LensPotentialFileName = concat(outroot,LensPotentialFileName)
        ScalarCovFileName =  Ini_Read_String_Default('scalar_covariance_output_file','scalCovCls.dat',.false.)
        if (ScalarCovFileName/='') then
            has_cl_2D_array = .true.
            ScalarCovFileName = concat(outroot,ScalarCovFileName)
        end if
    end if
    if (P%WantTensors) then
        TensorFileName =  trim(outroot) //Ini_Read_String('tensor_output_file')
        if (P%WantScalars)  then
            TotalFileName =  trim(outroot) //Ini_Read_String('total_output_file')
            LensedTotFileName = Ini_Read_String('lensed_total_output_file')
            if (LensedTotFileName/='') LensedTotFileName= trim(outroot) //trim(LensedTotFileName)
        end if
    end if
    if (P%WantVectors) then
        VectorFileName =  trim(outroot) //Ini_Read_String('vector_output_file')
    end if
    
!#SimDataAdd
    if(P%DoShePowFoC) then
        ShearFileName=trim(outroot)//Ini_Read_String('shear_output_file')
        GalPSFileName=trim(outroot)//Ini_Read_String('gal_output_file')
        ShearXGalPSFileName=trim(outroot)//Ini_Read_String('shear_x_gal_output_file')
        FracOfGalFileName=trim(outroot)//Ini_Read_String('frac_of_gal_file')
        RedShiftBinFileName=trim(outroot)//Ini_Read_String('redshift_bin_file')
        LMinValFileName=trim(outroot)//Ini_Read_String('lminval_file')
        LMaxValFileName=trim(outroot)//Ini_Read_String('lmaxval_file')
        ClFileName=trim(outroot)//Ini_Read_String('cl_output_file')
        ClRedshiftBinFileName=trim(outroot)//Ini_Read_String('cl_redshift_bin_file')
        ClMassBinFileName=trim(outroot)//Ini_Read_String('cl_mass_bin_file')
    end if
!#SimDataAdd    

#ifdef WRITE_FITS
    if (P%WantCls) then
        FITSfilename =  trim(outroot) //Ini_Read_String('FITS_filename',.true.)
        if (FITSfilename /='') then
            inquire(file=FITSfilename, exist=bad)
            if (bad) then
                open(unit=18,file=FITSfilename,status='old')
                close(18,status='delete')
            end if
        end if
    end if
#endif


    Ini_fail_on_not_found = .false.

    !optional parameters controlling the computation

    P%AccuratePolarization = Ini_Read_Logical('accurate_polarization',.true.)
    P%AccurateReionization = Ini_Read_Logical('accurate_reionization',.false.)
    P%AccurateBB = Ini_Read_Logical('accurate_BB',.false.)
    if (P%AccurateBB .and. P%WantCls .and. (P%Max_l < 3500 .or. &
        (P%NonLinear/=NonLinear_lens .and. P%NonLinear/=NonLinear_both) .or. P%Max_eta_k < 18000)) &
        write(*,*) 'WARNING: for accurate lensing BB you need high l_max_scalar, k_eta_max_scalar and non-linear lensing'

    P%DerivedParameters = Ini_Read_Logical('derived_parameters',.true.)

    version_check = Ini_Read_String('version_check')
    if (version_check == '') then
        !tag the output used parameters .ini file with the version of CAMB being used now
        call TNameValueList_Add(DefIni%ReadValues, 'version_check', version)
    else if (version_check /= version) then
        write(*,*) 'WARNING: version_check does not match this CAMB version'
    end if
    !Mess here to fix typo with backwards compatibility
    if (Ini_HasKey('do_late_rad_trunction')) then
        DoLateRadTruncation = Ini_Read_Logical('do_late_rad_trunction',.true.)
        if (Ini_HasKey('do_late_rad_truncation')) error stop 'check do_late_rad_xxxx'
    else
        DoLateRadTruncation = Ini_Read_Logical('do_late_rad_truncation',.true.)
    end if

    if (HighAccuracyDefault) then
        DoTensorNeutrinos = .true.
    else
        DoTensorNeutrinos = Ini_Read_Logical('do_tensor_neutrinos',DoTensorNeutrinos )
    end if
    FeedbackLevel = Ini_Read_Int('feedback_level',FeedbackLevel)

    output_file_headers = Ini_Read_Logical('output_file_headers',output_file_headers)

    P%MassiveNuMethod  = Ini_Read_Int('massive_nu_approx',Nu_best)

    ThreadNum      = Ini_Read_Int('number_of_threads',ThreadNum)
    use_spline_template = Ini_Read_Logical('use_spline_template',use_spline_template)

    if (do_bispectrum) then
        lSampleBoost   = 50
    else
        lSampleBoost   = Ini_Read_Double('l_sample_boost',lSampleBoost)
    end if
    if (outroot /= '') then
        if (InputFile /= trim(outroot) //'params.ini') then
            call Ini_SaveReadValues(trim(outroot) //'params.ini',1)
        else
            write(*,*) 'Output _params.ini not created as would overwrite input'
        end if
    end if

    call Ini_Close

    if (.not. CAMB_ValidateParams(P)) error stop 'Stopped due to parameter error'

#ifdef RUNIDLE
    call SetIdle
#endif

    if (global_error_flag==0) call CAMB_GetResults(P)
    if (global_error_flag/=0) then
        write(*,*) 'Error result '//trim(global_error_message)
        error stop
    endif
!#SimDataReplace
!Luci
!    if (P%PK_WantTransfer) then
     if (P%PK_WantTransfer .and. .not. P%DoShePowFoC) then
!SimDataReplace
        call Transfer_SaveToFiles(MT,TransferFileNames)
        call Transfer_SaveMatterPower(MT,MatterPowerFileNames)
        call Transfer_output_sig8(MT)
    end if
    
!#SimDataAdd
    if (P%OutputSimDataFiles .and. P%DoCMB) then
        call OutputSimDataCMB(P)
    end if
!#SimDataAdd

    if (P%WantCls) then
        call output_cl_files(ScalarFileName, ScalarCovFileName, TensorFileName, TotalFileName, &
            LensedFileName, LensedTotFilename, output_factor)

        call output_lens_pot_files(LensPotentialFileName, output_factor)

        if (P%WantVectors) then
            call output_veccl_files(VectorFileName, output_factor)
        end if

#ifdef WRITE_FITS
        if (FITSfilename /= '') call WriteFitsCls(FITSfilename, CP%Max_l)
#endif
    end if

    call CAMB_cleanup
    stop

100 stop 'Must give num_massive number of integer physical neutrinos for each eigenstate'
    end program driver


#ifdef RUNIDLE
    !If in Windows and want to run with low priorty so can multitask
    subroutine SetIdle
    USE DFWIN
    Integer dwPriority
    Integer CheckPriority

    dwPriority = 64 ! idle priority
    CheckPriority = SetPriorityClass(GetCurrentProcess(), dwPriority)

    end subroutine SetIdle
#endif
