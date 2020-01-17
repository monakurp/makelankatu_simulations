!> @file chem_emissions_mod.f90
!--------------------------------------------------------------------------------!
! This file is part of PALM model system.
!
! PALM is free software: you can redistribute it and/or modify it under the
! terms of the GNU General Public License as published by the Free Software
! Foundation, either version 3 of the License, or (at your option) any later
! version.
!
! PALM is distributed in the hope that it will be useful, but WITHOUT ANY
! WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
! A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License along with
! PALM. If not, see <http://www.gnu.org/licenses/>.
!
! Copyright 2018-2019 Leibniz Universitaet Hannover
! Copyright 2018-2019 Freie Universitaet Berlin
! Copyright 2018-2019 Karlsruhe Institute of Technology
!--------------------------------------------------------------------------------!
!
! Current revisions:
! ------------------
! 
! 
! Former revisions:
! -----------------
! $Id: chem_emissions_mod.f90 4356 2019-12-20 17:09:33Z suehring $
! Minor formatting adjustment
! 
! 4242 2019-09-27 12:59:10Z suehring
! Adjust index_hh access to new definition accompanied with new 
! palm_date_time_mod. Note, this is just a preleminary fix. (E Chan)
! 
! 4230 2019-09-11 13:58:14Z suehring
! Bugfix, consider that time_since_reference_point can be also negative when 
! time indices are determined.
! 
! 4227 2019-09-10 18:04:34Z gronemeier
! implement new palm_date_time_mod
! 
! 4223 2019-09-10 09:20:47Z gronemeier
! Unused routine chem_emissions_check_parameters commented out due to uninitialized content
! 
! 4182 2019-08-22 15:20:23Z scharf
! Corrected "Former revisions" section
! 
! 4157 2019-08-14 09:19:12Z suehring
! Replace global arrays also in mode_emis branch
! 
! 4154 2019-08-13 13:35:59Z suehring
! Replace global arrays for emissions by local ones.
! 
! 4144 2019-08-06 09:11:47Z raasch
! relational operators .EQ., .NE., etc. replaced by ==, /=, etc.
! 
! 4055 2019-06-27 09:47:29Z suehring
! - replaced local thermo. constants w/ module definitions in
!   basic_constants_and_equations_mod (rgas_univ, p_0, r_d_cp) 
! - initialize array emis_distribution immediately following allocation
! - lots of minor formatting changes based on review sesson in 20190325 
!   (E.C. Chan)
! 
! 3968 2019-05-13 11:04:01Z suehring
! - in subroutine chem_emissions_match replace all decision structures relating to
!   mode_emis to emiss_lod
! - in subroutine chem_check_parameters replace emt%nspec with emt%n_emiss_species
! - spring cleaning (E.C. Chan)
! 
! 3885 2019-04-11 11:29:34Z kanani
! Changes related to global restructuring of location messages and introduction 
! of additional debug messages
! 
! 3831 2019-03-28 09:11:22Z forkel
! added nvar to USE chem_gasphase_mod (chem_modules will not include nvar anymore) 
! 
! 3788 2019-03-07 11:40:09Z banzhafs
! Removed unused variables from chem_emissions_mod
! 
! 3772 2019-02-28 15:51:57Z suehring 
! - In case of parametrized emissions, assure that emissions are only on natural
!   surfaces (i.e. streets) and not on urban surfaces. 
! - some unnecessary if clauses removed
!
! 3685 2019 -01-21 01:02:11Z knoop
! Some interface calls moved to module_interface + cleanup
! 3286 2018-09-28 07:41:39Z forkel
!
! Authors:
! --------
! @author Emmanuele Russo (FU-Berlin)
! @author Sabine Banzhaf  (FU-Berlin)
! @author Martijn Schaap  (FU-Berlin, TNO Utrecht)
! 
! Description:
! ------------
!>  MODULE for reading-in Chemistry Emissions data
!>
!> @todo Rename nspec to n_emis to avoid inteferece with nspec from chem_gasphase_mod
!> @todo Check_parameters routine should be developed: for now it includes just one condition
!> @todo Use of Restart files not contempled at the moment
!> @todo revise indices of files read from the netcdf: preproc_emission_data and expert_emission_data
!> @todo for now emission data may be passed on a singular vertical level: need to be more flexible
!> @todo fill/activate restart option in chem_emissions_init
!> @todo discuss dt_emis
!> @note <Enter notes on the module>
!> @bug  <Enter known bugs here>
!------------------------------------------------------------------------------!

 MODULE chem_emissions_mod

    USE arrays_3d,                                                          &
        ONLY:  rho_air

    USE basic_constants_and_equations_mod,                                  &
        ONLY:  rgas_univ, p_0, rd_d_cp

    USE control_parameters,                                                 &
        ONLY:  debug_output,                                                &
               end_time, message_string, initializing_actions,              &
               intermediate_timestep_count, dt_3d
 
    USE indices

    USE kinds

#if defined( __netcdf )
    USE netcdf
#endif

    USE netcdf_data_input_mod,                                               &
        ONLY: chem_emis_att_type, chem_emis_val_type

    USE chem_gasphase_mod,                                                   &
        ONLY: nvar, spc_names
 
    USE chem_modules

    USE statistics,                                                          &   
        ONLY:  weight_pres

    
    IMPLICIT NONE

!
!-- Declare all global variables within the module 
    
    CHARACTER (LEN=80) ::  filename_emis             !< Variable for the name of the netcdf input file

    INTEGER(iwp) ::  dt_emis                         !< Time Step Emissions 
    INTEGER(iwp) ::  i                               !< index 1st selected dimension (some dims are not spatial)
    INTEGER(iwp) ::  j                               !< index 2nd selected dimension 
    INTEGER(iwp) ::  i_start                         !< Index to start read variable from netcdf along one dims
    INTEGER(iwp) ::  i_end                           !< Index to end read variable from netcdf in one dims
    INTEGER(iwp) ::  j_start                         !< Index to start read variable from netcdf in additional dims
    INTEGER(iwp) ::  j_end                           !< Index to end read variable from netcdf in additional dims
    INTEGER(iwp) ::  len_index                       !< length of index (used for several indices)
    INTEGER(iwp) ::  len_index_pm                    !< length of PMs index
    INTEGER(iwp) ::  len_index_voc                   !< length of voc index
    INTEGER(iwp) ::  z_start                         !< Index to start read variable from netcdf in additional dims
    INTEGER(iwp) ::  z_end                           !< Index to end read variable from netcdf in additional dims

    REAL(wp) ::  conversion_factor                   !< Units Conversion Factor

    SAVE

! !
! !-- Checks Input parameters
!     INTERFACE chem_emissions_check_parameters
!        MODULE PROCEDURE chem_emissions_check_parameters
!     END INTERFACE chem_emissions_check_parameters
!
!-- Matching Emissions actions  
    INTERFACE chem_emissions_match
       MODULE PROCEDURE chem_emissions_match
    END INTERFACE chem_emissions_match
!
!-- Initialization actions  
    INTERFACE chem_emissions_init
       MODULE PROCEDURE chem_emissions_init
    END INTERFACE chem_emissions_init
!
!-- Setup of Emissions
    INTERFACE chem_emissions_setup
       MODULE PROCEDURE chem_emissions_setup
    END INTERFACE chem_emissions_setup
    
    PUBLIC chem_emissions_init, chem_emissions_match, chem_emissions_setup
!
!-- Public Variables
    PUBLIC conversion_factor, len_index, len_index_pm, len_index_voc

 CONTAINS

! !------------------------------------------------------------------------------!
! ! Description:
! ! ------------
! !> Routine for checking input parameters
! !------------------------------------------------------------------------------!
!  SUBROUTINE chem_emissions_check_parameters
!
!     IMPLICIT NONE
!
!     TYPE(chem_emis_att_type) ::  emt
!
! !
! !-- Check if species count matches the number of names
! !-- passed for the chemiscal species
!
!     IF  ( SIZE(emt%species_name) /= emt%n_emiss_species  )  THEN
! !    IF  ( SIZE(emt%species_name) /= emt%nspec  )  THEN
!
!        message_string = 'Numbers of input emission species names and number of species'     //      &
!                          'for which emission values are given do not match'
!        CALL message( 'chem_emissions_check_parameters', 'CM0437', 2, 2, 0, 6, 0 )
!
!     ENDIF
!
!  END SUBROUTINE chem_emissions_check_parameters


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Matching the chemical species indices. The routine checks what are the
!> indices of the emission input species and the corresponding ones of the
!> model species. The routine gives as output a vector containing the number 
!> of common species: it is important to note that while the model species
!> are distinct, their values could be given to a single species in input.
!> For example, in the case of NO2 and NO, values may be passed in input as
!> NOX values.
!------------------------------------------------------------------------------!

SUBROUTINE chem_emissions_match( emt_att,len_index )   

    INTEGER(iwp)  ::  ind_inp                    !< Parameters for cycling through chemical input species 
    INTEGER(iwp)  ::  ind_mod                    !< Parameters for cycling through chemical model species 
    INTEGER(iwp)  ::  ind_voc                    !< Indices to check whether a split for voc should be done
    INTEGER(iwp)  ::  ispec                      !< index for cycle over effective number of emission species
    INTEGER(iwp)  ::  nspec_emis_inp             !< Variable where to store # of emission species in input

    INTEGER(iwp), INTENT(INOUT)  ::  len_index   !< number of common species between input dataset & model species

    TYPE(chem_emis_att_type), INTENT(INOUT) ::  emt_att     !< Chemistry Emission Array (decl. netcdf_data_input.f90)


    IF  ( debug_output )  CALL debug_message( 'chem_emissions_match', 'start' )

!
!-- Number of input emission species

    nspec_emis_inp = emt_att%n_emiss_species
!    nspec_emis_inp=emt_att%nspec 

!
!-- Check the emission LOD: 0 (PARAMETERIZED), 1 (DEFAULT), 2 (PRE-PROCESSED)
!
    SELECT CASE (emiss_lod)

!
!-- LOD 0 (PARAMETERIZED mode)

       CASE (0)

          len_index = 0

! number of species and number of matched species can be different
! but call is only made if both are greater than zero

          IF  ( nvar > 0  .AND.  nspec_emis_inp > 0 )  THEN

!
!-- Cycle over model species

             DO  ind_mod = 1, nvar
                ind_inp = 1
                DO  WHILE ( TRIM( surface_csflux_name(ind_inp) ) /= 'novalue' )    !< 'novalue' is the default  
                   IF  ( TRIM( surface_csflux_name(ind_inp) ) == TRIM( spc_names(ind_mod) ) )  THEN
                      len_index = len_index + 1
                   ENDIF
                   ind_inp = ind_inp + 1
                ENDDO
             ENDDO

             IF  ( len_index > 0 )  THEN

!
!-- Allocation of Arrays of the matched species

                ALLOCATE ( match_spec_input(len_index) ) 
                ALLOCATE ( match_spec_model(len_index) )

!
!-- Pass species indices to declared arrays

                len_index = 0

                DO  ind_mod = 1, nvar  
                   ind_inp = 1
                   DO  WHILE  ( TRIM( surface_csflux_name(ind_inp) ) /= 'novalue' )
                      IF  ( TRIM(surface_csflux_name(ind_inp)) ==          &
                            TRIM(spc_names(ind_mod))            )  THEN
                         len_index = len_index + 1
                         match_spec_input(len_index) = ind_inp
                         match_spec_model(len_index) = ind_mod
                      ENDIF
                   ind_inp = ind_inp + 1
                   END DO
                END DO

!
!-- Check

                DO  ispec = 1, len_index

                   IF  ( emiss_factor_main(match_spec_input(ispec) ) < 0    .AND.     &
                         emiss_factor_side(match_spec_input(ispec) ) < 0 )  THEN

                      message_string = 'PARAMETERIZED emissions mode selected:'             //          &
                                       ' EMISSIONS POSSIBLE ONLY ON STREET SURFACES'        //          &
                                       ' but values of scaling factors for street types'    //          &
                                       ' emiss_factor_main AND emiss_factor_side'           //          &
                                       ' not provided for each of the emissions species'    //          &
                                       ' or not provided at all: PLEASE set a finite value' //          &
                                       ' for these parameters in the chemistry namelist'          
                      CALL message( 'chem_emissions_matching', 'CM0442', 2, 2, 0, 6, 0 )
 
                   ENDIF

                END DO


             ELSE
                
                message_string = 'Non of given Emission Species'           //           &
                                 ' matches'                                //           &
                                 ' model chemical species'                 //           &
                                 ' Emission routine is not called'         
                CALL message( 'chem_emissions_matching', 'CM0443', 0, 0, 0, 6, 0 )
 
             ENDIF

          ELSE
     
             message_string = 'Array of Emission species not allocated: '             //          &
                              ' Either no emission species are provided as input or'  //          &
                              ' no chemical species are used by PALM.'                //          &
                              ' Emission routine is not called'                                   
             CALL message( 'chem_emissions_matching', 'CM0444', 0, 2, 0, 6, 0 ) 
 
          ENDIF

!
!-- LOD 1 (DEFAULT mode)

       CASE (1)

          len_index = 0          ! total number of species (to be accumulated)  
          len_index_voc = 0      ! total number of VOCs (to be accumulated)
          len_index_pm = 3       ! total number of PMs: PM1, PM2.5, PM10.

!
!-- number of model species and input species could be different
!-- but process this only when both are non-zero
 
          IF  ( nvar > 0  .AND.  nspec_emis_inp > 0 )  THEN

!
!-- Cycle over model species
             DO  ind_mod = 1, nvar

!
!-- Cycle over input species

                DO  ind_inp = 1, nspec_emis_inp

!
!-- Check for VOC Species

                   IF  ( TRIM( emt_att%species_name(ind_inp) ) == "VOC" )  THEN
                      DO  ind_voc= 1, emt_att%nvoc
                        
                         IF  ( TRIM( emt_att%voc_name(ind_voc) ) == TRIM( spc_names(ind_mod) ) )  THEN
                            len_index = len_index + 1
                            len_index_voc = len_index_voc + 1
                         ENDIF
                         
                      END DO
                   ENDIF

!
!-- PMs: There is one input species name for all PM 
!-- This variable has 3 dimensions, one for PM1, PM2.5 and PM10

                   IF  ( TRIM( emt_att%species_name(ind_inp) ) == "PM" )  THEN

                      IF      ( TRIM( spc_names(ind_mod) ) == "PM1" )   THEN
                         len_index = len_index + 1
                      ELSEIF  ( TRIM( spc_names(ind_mod) ) == "PM25" )  THEN
                         len_index = len_index + 1
                      ELSEIF  ( TRIM( spc_names(ind_mod) ) == "PM10" )  THEN
                         len_index = len_index + 1
                      ENDIF

                   ENDIF

!
!-- NOX: NO2 and NO

                   IF  ( TRIM( emt_att%species_name(ind_inp) ) == "NOX" )  THEN 

                      IF     ( TRIM( spc_names(ind_mod) ) == "NO"  )  THEN
                         len_index = len_index + 1
                      ELSEIF  ( TRIM( spc_names(ind_mod) ) == "NO2" )  THEN
                         len_index = len_index + 1
                      ENDIF

                   ENDIF

!
!-- SOX: SO2 and SO4

                   IF  ( TRIM( emt_att%species_name(ind_inp) ) == "SOX" )  THEN

                      IF      ( TRIM( spc_names(ind_mod) ) == "SO2" )  THEN
                         len_index = len_index + 1
                      ELSEIF  ( TRIM( spc_names(ind_mod) ) == "SO4" )  THEN
                         len_index = len_index + 1
                      ENDIF

                   ENDIF

!
!-- Other Species

                   IF  ( TRIM( emt_att%species_name(ind_inp) ) ==             &
                        TRIM( spc_names(ind_mod) ) )  THEN
                      len_index = len_index + 1
                   ENDIF

                END DO  ! ind_inp ...

             END DO  ! ind_mod ...


!
!-- Allocate arrays

             IF  ( len_index > 0 )  THEN

                ALLOCATE ( match_spec_input(len_index) )
                ALLOCATE ( match_spec_model(len_index) )

                IF  ( len_index_voc > 0 )  THEN

!
!-- Contains indices of the VOC model species

                   ALLOCATE( match_spec_voc_model(len_index_voc) )

!
!-- Contains the indices of different values of VOC composition
!-- of input variable VOC_composition

                   ALLOCATE( match_spec_voc_input(len_index_voc) )

                ENDIF

!
!-- Pass the species indices to declared arrays

                len_index = 0
                len_index_voc = 0
                
                DO  ind_mod = 1, nvar
                   DO  ind_inp = 1, nspec_emis_inp 

!
!-- VOCs

                      IF  ( TRIM( emt_att%species_name(ind_inp) ) == "VOC"  .AND.        &
                            ALLOCATED (match_spec_voc_input) )  THEN

                         DO  ind_voc = 1, emt_att%nvoc

                            IF  ( TRIM( emt_att%voc_name(ind_voc) ) ==                  &
                                 TRIM( spc_names(ind_mod) ) )  THEN

                               len_index     = len_index + 1
                               len_index_voc = len_index_voc + 1
                        
                               match_spec_input(len_index) = ind_inp
                               match_spec_model(len_index) = ind_mod

                               match_spec_voc_input(len_index_voc) = ind_voc
                               match_spec_voc_model(len_index_voc) = ind_mod

                            ENDIF

                         END DO

                      ENDIF

!
!-- PMs

                      IF  ( TRIM( emt_att%species_name(ind_inp) ) == "PM" )  THEN

                         IF      ( TRIM( spc_names(ind_mod) ) == "PM1"  )  THEN
                            len_index = len_index + 1
                            match_spec_input(len_index) = ind_inp
                            match_spec_model(len_index) = ind_mod
                         ELSEIF  ( TRIM( spc_names(ind_mod) ) == "PM25" )  THEN
                            len_index = len_index + 1
                            match_spec_input(len_index) = ind_inp
                            match_spec_model(len_index) = ind_mod
                         ELSEIF  ( TRIM( spc_names(ind_mod) ) == "PM10" )  THEN
                            len_index = len_index + 1
                            match_spec_input(len_index) = ind_inp
                            match_spec_model(len_index) = ind_mod
                         ENDIF

                      ENDIF

!
!-- NOX 
                      IF  ( TRIM( emt_att%species_name(ind_inp) ) == "NOX" )  THEN

                         IF      ( TRIM( spc_names(ind_mod) ) == "NO"  )  THEN
                            len_index = len_index + 1

                            match_spec_input(len_index) = ind_inp
                            match_spec_model(len_index) = ind_mod

                         ELSEIF  ( TRIM( spc_names(ind_mod) ) == "NO2" )  THEN
                            len_index = len_index + 1

                            match_spec_input(len_index) = ind_inp
                            match_spec_model(len_index) = ind_mod
  
                         ENDIF

                      ENDIF


!
!-- SOX 

                      IF  ( TRIM( emt_att%species_name(ind_inp) ) == "SOX" ) THEN

                         IF  ( TRIM( spc_names(ind_mod) ) == "SO2" )  THEN
                            len_index = len_index + 1
                            match_spec_input(len_index) = ind_inp
                            match_spec_model(len_index) = ind_mod
                         ELSEIF  ( TRIM( spc_names(ind_mod) ) == "SO4" )  THEN
                            len_index = len_index + 1
                            match_spec_input(len_index) = ind_inp
                            match_spec_model(len_index) = ind_mod
                         ENDIF

                      ENDIF

!
!-- Other Species

                      IF  ( TRIM( emt_att%species_name(ind_inp) ) == TRIM( spc_names(ind_mod) ) )  THEN
                         len_index = len_index + 1
                         match_spec_input(len_index) = ind_inp
                         match_spec_model(len_index) = ind_mod
                      ENDIF

                   END DO  ! inp_ind

                END DO     ! inp_mod

!
!-- Error reporting (no matching)

             ELSE

                message_string = 'None of given Emission Species matches'           //   &
                                 ' model chemical species'                          //   &
                                 ' Emission routine is not called'         
                CALL message( 'chem_emissions_matching', 'CM0440', 0, 0, 0, 6, 0 ) 

             ENDIF

!
!-- Error reporting (no species)

          ELSE

             message_string = 'Array of Emission species not allocated: '             // &
                              ' Either no emission species are provided as input or'  // &
                              ' no chemical species are used by PALM:'                // &
                              ' Emission routine is not called'                                   
             CALL message( 'chem_emissions_matching', 'CM0441', 0, 2, 0, 6, 0 ) 
 
          ENDIF

!
!-- LOD 2 (PRE-PROCESSED mode)

       CASE (2)

          len_index = 0
          len_index_voc = 0

          IF  ( nvar > 0  .AND.  nspec_emis_inp > 0 )  THEN
!
!-- Cycle over model species
             DO  ind_mod = 1, nvar 

!
!-- Cycle over input species  
                DO  ind_inp = 1, nspec_emis_inp

!
!-- Check for VOC Species 

                   IF  ( TRIM( emt_att%species_name(ind_inp) ) == "VOC" )  THEN        
                      DO  ind_voc = 1, emt_att%nvoc
                         IF  ( TRIM( emt_att%voc_name(ind_voc) ) == TRIM( spc_names(ind_mod) ) )  THEN
                            len_index     = len_index + 1
                            len_index_voc = len_index_voc + 1
                         ENDIF
                      END DO
                   ENDIF

!
!-- Other Species

                   IF  ( TRIM(emt_att%species_name(ind_inp)) == TRIM(spc_names(ind_mod)) )  THEN
                      len_index = len_index + 1
                   ENDIF
                ENDDO
             ENDDO

!
!-- Allocate array for storing the indices of the matched species

             IF  ( len_index > 0 )  THEN
 
                ALLOCATE ( match_spec_input(len_index) ) 
 
                ALLOCATE ( match_spec_model(len_index) )

                IF  ( len_index_voc > 0 )  THEN
!
!-- contains indices of the VOC model species
                   ALLOCATE( match_spec_voc_model(len_index_voc) )
!
!-- contains the indices of different values of VOC composition of input variable VOC_composition
                   ALLOCATE( match_spec_voc_input(len_index_voc) )

                ENDIF

!
!-- pass the species indices to declared arrays

                len_index = 0

!
!-- Cycle over model species

                DO  ind_mod = 1, nvar
 
!
!-- Cycle over Input species  

                   DO  ind_inp = 1, nspec_emis_inp

!
!-- VOCs

                      IF  ( TRIM(emt_att%species_name(ind_inp) ) == "VOC"  .AND.     &
                           ALLOCATED(match_spec_voc_input) )  THEN
                         
                         DO  ind_voc= 1, emt_att%nvoc
                            IF  ( TRIM( emt_att%voc_name(ind_voc) ) == TRIM( spc_names(ind_mod) ) )  THEN
                               len_index = len_index + 1
                               len_index_voc = len_index_voc + 1
                        
                               match_spec_input(len_index) = ind_inp
                               match_spec_model(len_index) = ind_mod

                               match_spec_voc_input(len_index_voc) = ind_voc
                               match_spec_voc_model(len_index_voc) = ind_mod                         
                            ENDIF
                         END DO
                      ENDIF

!
!-- Other Species

                      IF  ( TRIM( emt_att%species_name(ind_inp) ) == TRIM( spc_names(ind_mod) ) )  THEN
                         len_index = len_index + 1
                         match_spec_input(len_index) = ind_inp
                         match_spec_model(len_index) = ind_mod
                      ENDIF

                   END DO  ! ind_inp
                END DO     ! ind_mod

             ELSE  ! if len_index_voc <= 0

!
!-- in case there are no species matching (just informational message)

                message_string = 'Non of given emission species'            //         &
                                 ' matches'                                //          &
                                 ' model chemical species:'                //          &
                                 ' Emission routine is not called'  
                CALL message( 'chem_emissions_matching', 'CM0438', 0, 0, 0, 6, 0 )
             ENDIF

!
!-- Error check (no matching)
 
          ELSE

!
!-- either spc_names is zero or nspec_emis_inp is not allocated
             message_string = 'Array of Emission species not allocated:'              // &
                              ' Either no emission species are provided as input or'  // &
                              ' no chemical species are used by PALM:'                // &
                              ' Emission routine is not called'                  
             CALL message( 'chem_emissions_matching', 'CM0439', 0, 2, 0, 6, 0 ) 

          ENDIF  

!
!-- If emission module is switched on but mode_emis is not specified or it is given the wrong name

!
!-- Error check (no species)

       CASE DEFAULT

          message_string = 'Emission Module switched ON, but'                           //         &
                           ' either no emission mode specified or incorrectly given :'  //         &
                           ' please, pass the correct value to the namelist parameter "mode_emis"'                  
          CALL message( 'chem_emissions_matching', 'CM0445', 2, 2, 0, 6, 0 )              

       END SELECT

       IF  ( debug_output )  CALL debug_message( 'chem_emissions_match', 'end' )

 END SUBROUTINE chem_emissions_match

 
!------------------------------------------------------------------------------!
! Description:
! ------------
!> Initialization:
!> Netcdf reading, arrays allocation and first calculation of cssws
!> fluxes at timestep 0 
!------------------------------------------------------------------------------!

 SUBROUTINE chem_emissions_init

    USE netcdf_data_input_mod,                                              &
        ONLY:  chem_emis, chem_emis_att
    
    IMPLICIT NONE
 
    INTEGER(iwp) :: ispec                        !< running index

!    
!-- Actions for initial runs 
!  IF  (  TRIM( initializing_actions ) /= 'read_restart_data' )  THEN
!--    ... 
!    
!
!-- Actions for restart runs
!  ELSE
!--    ...
!
!  ENDIF


    IF  ( debug_output )  CALL debug_message( 'chem_emissions_init', 'start' )

!
!-- Matching

    CALL  chem_emissions_match( chem_emis_att, n_matched_vars )
 
    IF  ( n_matched_vars == 0 )  THEN
  
       emission_output_required = .FALSE.

    ELSE

       emission_output_required = .TRUE.


!
!-- Set molecule masses  (in kg/mol)

       ALLOCATE( chem_emis_att%xm(n_matched_vars) )

       DO  ispec = 1, n_matched_vars
          SELECT CASE ( TRIM( spc_names(match_spec_model(ispec)) ) )
             CASE ( 'SO2'  );  chem_emis_att%xm(ispec) = xm_S + xm_O * 2
             CASE ( 'SO4'  );  chem_emis_att%xm(ispec) = xm_S + xm_O * 4
             CASE ( 'NO'   );  chem_emis_att%xm(ispec) = xm_N + xm_O
             CASE ( 'NO2'  );  chem_emis_att%xm(ispec) = xm_N + xm_O * 2
             CASE ( 'NH3'  );  chem_emis_att%xm(ispec) = xm_N + xm_H * 3
             CASE ( 'CO'   );  chem_emis_att%xm(ispec) = xm_C + xm_O
             CASE ( 'CO2'  );  chem_emis_att%xm(ispec) = xm_C + xm_O * 2
             CASE ( 'CH4'  );  chem_emis_att%xm(ispec) = xm_C + xm_H * 4
             CASE ( 'HNO3' );  chem_emis_att%xm(ispec) = xm_H + xm_N + xm_O*3
!***************************************************************************************************
! MONA ADDED:
             CASE ( 'OCSV' );  chem_emis_att%xm(ispec) = 150.0_wp
             CASE ( 'RH' );    chem_emis_att%xm(ispec) = 150.0_wp
             CASE ( 'H2SO4' ); chem_emis_att%xm(ispec) = xm_H * 2 + xm_S + xm_O * 4
!***************************************************************************************************
             CASE DEFAULT
                chem_emis_att%xm(ispec) = 1.0_wp
          END SELECT
       ENDDO

    
!
!-- Get emissions for the first time step base on LOD (if defined)
!-- or emission mode (if no LOD defined)

!
!-- NOTE - I could use a combined if ( lod = xxx .or. mode = 'XXX' )
!--        type of decision structure but I think it is much better
!--        to implement it this way (i.e., conditional on lod if it
!--        is defined, and mode if not) as we can easily take out
!--        the case structure for mode_emis later on.

       IF   ( emiss_lod < 0 )  THEN  !-- no LOD defined (not likely)

          SELECT CASE ( TRIM( mode_emis ) )   

             CASE ( 'PARAMETERIZED' )     ! LOD 0

                IF  (  .NOT.  ALLOCATED( emis_distribution) )  THEN
                   ALLOCATE( emis_distribution(1,nys:nyn,nxl:nxr,n_matched_vars) )
                ENDIF

                CALL chem_emissions_setup( chem_emis_att, chem_emis, n_matched_vars)

             CASE ( 'DEFAULT' )           ! LOD 1

                IF  (  .NOT.  ALLOCATED( emis_distribution) )  THEN
                   ALLOCATE( emis_distribution(1,nys:nyn,nxl:nxr,n_matched_vars) )
                ENDIF

                CALL chem_emissions_setup( chem_emis_att, chem_emis, n_matched_vars )

             CASE ( 'PRE-PROCESSED' )     ! LOD 2

                IF  (  .NOT.  ALLOCATED( emis_distribution) )  THEN
!
!--                Note, at the moment emissions are considered only by surface fluxes
!--                rather than by volume sources. Therefore, no vertical dimension is 
!--                required and is thus allocated with 1. Later when volume sources
!--                are considered, the vertical dimension will increase.                 
                   !ALLOCATE( emis_distribution(nzb:nzt+1,nys:nyn,nxl:nxr,n_matched_vars) )
                   ALLOCATE( emis_distribution(1,nys:nyn,nxl:nxr,n_matched_vars) )
                ENDIF
 
                CALL chem_emissions_setup( chem_emis_att, chem_emis, n_matched_vars )

          END SELECT

       ELSE  ! if LOD is defined

          SELECT CASE ( emiss_lod )

             CASE ( 0 )     ! parameterized mode

                IF  (  .NOT.  ALLOCATED( emis_distribution) )  THEN
                   ALLOCATE( emis_distribution(1,nys:nyn,nxl:nxr,n_matched_vars) )
                ENDIF

                CALL chem_emissions_setup( chem_emis_att, chem_emis, n_matched_vars)

             CASE ( 1 )     ! default mode

                IF  (  .NOT.  ALLOCATED( emis_distribution) )  THEN
                   ALLOCATE( emis_distribution(1,nys:nyn,nxl:nxr,n_matched_vars) )
                ENDIF

                CALL chem_emissions_setup( chem_emis_att, chem_emis, n_matched_vars )

             CASE ( 2 )     ! pre-processed mode

                IF  (  .NOT.  ALLOCATED( emis_distribution) )  THEN
                   ALLOCATE( emis_distribution(nzb:nzt+1,nys:nyn,nxl:nxr,n_matched_vars) )
                ENDIF
 
                CALL chem_emissions_setup( chem_emis_att, chem_emis, n_matched_vars )

          END SELECT

       ENDIF

!
! -- initialize

       emis_distribution = 0.0_wp

    ENDIF

    IF  ( debug_output )  CALL debug_message( 'chem_emissions_init', 'end' )

 END SUBROUTINE chem_emissions_init



!------------------------------------------------------------------------------!
! Description:
! ------------
!> Routine for Update of Emission values at each timestep.
!>
!> @todo Clarify the correct usage of index_dd, index_hh and index_mm. Consider
!>       renaming of these variables.
!> @todo Clarify time used in emis_lod=2 mode. ATM, the used time seems strange.
!-------------------------------------------------------------------------------!

 SUBROUTINE chem_emissions_setup( emt_att, emt, n_matched_vars )
 
   USE surface_mod,                                               &
       ONLY:  surf_def_h, surf_lsm_h, surf_usm_h

   USE netcdf_data_input_mod,                                     &
       ONLY:  street_type_f

   USE arrays_3d,                                                 &        
       ONLY: hyp, pt 

    USE control_parameters, &
        ONLY: time_since_reference_point

    USE palm_date_time_mod, &
        ONLY: days_per_week, get_date_time, hours_per_day, months_per_year, seconds_per_day
   
    IMPLICIT NONE
  

    TYPE(chem_emis_att_type), INTENT(INOUT) ::  emt_att                         !< variable to store emission information 

    TYPE(chem_emis_val_type), INTENT(INOUT), ALLOCATABLE, DIMENSION(:) ::  emt  !< variable to store emission input values, 
                                                                                !< depending on the considered species

    INTEGER,INTENT(IN) ::  n_matched_vars                                       !< Output of matching routine with number
                                                                                !< of matched species

    INTEGER(iwp) ::  i                                                          !< running index for grid in x-direction
    INTEGER(iwp) ::  i_pm_comp                                                  !< index for number of PM components
    INTEGER(iwp) ::  icat                                                       !< Index for number of categories
    INTEGER(iwp) ::  ispec                                                      !< index for number of species
    INTEGER(iwp) ::  ivoc                                                       !< Index for number of VOCs
    INTEGER(iwp) ::  j                                                          !< running index for grid in y-direction
    INTEGER(iwp) ::  k                                                          !< running index for grid in z-direction
    INTEGER(iwp) ::  m                                                          !< running index for horizontal surfaces

    INTEGER(iwp) ::  day_of_month                                               !< day of the month
    INTEGER(iwp) ::  day_of_week                                                !< day of the week
    INTEGER(iwp) ::  day_of_year                                                !< day of the year
    INTEGER(iwp) ::  days_since_reference_point                                 !< days since reference point
    INTEGER(iwp) ::  hour_of_day                                                !< hour of the day
    INTEGER(iwp) ::  month_of_year                                              !< month of the year
    INTEGER(iwp) ::  index_dd                                                   !< index day
    INTEGER(iwp) ::  index_hh                                                   !< index hour
    INTEGER(iwp) ::  index_mm                                                   !< index month

    REAL(wp) ::  time_utc_init                                                  !< second of day of initial date

    !
    !-- CONVERSION FACTORS: TIME  
    REAL(wp), PARAMETER ::  hour_per_year =  8760.0_wp  !< number of hours in a year of 365 days  
    REAL(wp), PARAMETER ::  s_per_hour    =  3600.0_wp  !< number of sec per hour (s)/(hour)    
    REAL(wp), PARAMETER ::  s_per_day     = 86400.0_wp  !< number of sec per day (s)/(day)  

    REAL(wp), PARAMETER ::  day_to_s      = 1.0_wp/s_per_day                   !< conversion day   -> sec
    REAL(wp), PARAMETER ::  hour_to_s     = 1.0_wp/s_per_hour                  !< conversion hours -> sec
    REAL(wp), PARAMETER ::  year_to_s     = 1.0_wp/(s_per_hour*hour_per_year)  !< conversion year  -> sec
    !
    !-- CONVERSION FACTORS: MASS
    REAL(wp), PARAMETER ::  g_to_kg       = 1.0E-03_wp     !< Conversion from g to kg (kg/g) 
    REAL(wp), PARAMETER ::  miug_to_kg    = 1.0E-09_wp     !< Conversion from g to kg (kg/g)
    REAL(wp), PARAMETER ::  tons_to_kg    = 100.0_wp       !< Conversion from tons to kg (kg/tons)    
    !
    !-- CONVERSION FACTORS: PPM
    REAL(wp), PARAMETER ::  ratio2ppm     = 1.0E+06_wp 
 
    REAL(wp), DIMENSION(24)                        ::  par_emis_time_factor    !< time factors for the parameterized mode:
                                                                               !< fixed houlry profile for example day
    REAL(wp), DIMENSION(nzb:nzt+1,nys:nyn,nxl:nxr) ::  conv_to_ratio           !< factor used for converting input 
                                                                               !< to concentration ratio
    REAL(wp), DIMENSION(nzb:nzt+1,nys:nyn,nxl:nxr) ::  tmp_temp                !< temporary variable for abs. temperature

    REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  delta_emis   !< incremental emission factor                     
    REAL(wp), DIMENSION(:),   ALLOCATABLE ::  time_factor  !< factor for time scaling of emissions
    REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  emis         !< emission factor

    IF  ( emission_output_required )  THEN

!
!-- Set emis_dt to be used - since chemistry ODEs can be stiff, the option
!-- to solve them at every RK substep is present to help improve stability
!-- should the need arises
 
       IF  ( call_chem_at_all_substeps )  THEN

          dt_emis = dt_3d * weight_pres(intermediate_timestep_count)

       ELSE

          dt_emis = dt_3d

       ENDIF

 !
 !-- Conversion of units to the ones employed in PALM 
 !-- In PARAMETERIZED mode no conversion is performed: in this case input units are fixed

        IF  ( TRIM( mode_emis ) == "DEFAULT"  .OR.  TRIM( mode_emis ) == "PRE-PROCESSED" )  THEN

          SELECT CASE ( TRIM( emt_att%units ) )

             CASE ( 'kg/m2/s', 'KG/M2/S'    );     conversion_factor = 1.0_wp                   ! kg
             CASE ( 'kg/m2/hour', 'KG/M2/HOUR' );  conversion_factor = hour_to_s
             CASE ( 'kg/m2/day', 'KG/M2/DAY'  );   conversion_factor = day_to_s
             CASE ( 'kg/m2/year', 'KG/M2/YEAR' );  conversion_factor = year_to_s

             CASE ( 'ton/m2/s', 'TON/M2/S'    );     conversion_factor = tons_to_kg             ! tonnes
             CASE ( 'ton/m2/hour', 'TON/M2/HOUR' );  conversion_factor = tons_to_kg*hour_to_s
             CASE ( 'ton/m2/year', 'TON/M2/YEAR' );  conversion_factor = tons_to_kg*year_to_s

             CASE ( 'g/m2/s', 'G/M2/S'    );     conversion_factor = g_to_kg                    ! grams
             CASE ( 'g/m2/hour', 'G/M2/HOUR' );  conversion_factor = g_to_kg*hour_to_s
             CASE ( 'g/m2/year', 'G/M2/YEAR' );  conversion_factor = g_to_kg*year_to_s

             CASE ( 'micrograms/m2/s', 'MICROGRAMS/M2/S'    );     conversion_factor = miug_to_kg            ! ug
             CASE ( 'micrograms/m2/hour', 'MICROGRAMS/M2/HOUR' );  conversion_factor = miug_to_kg*hour_to_s
             CASE ( 'micrograms/m2/year', 'MICROGRAMS/M2/YEAR' );  conversion_factor = miug_to_kg*year_to_s

!
!-- Error check (need units)

             CASE DEFAULT  
                message_string = 'The Units of the provided emission input'                 // &
                                 ' are not the ones required by PALM-4U: please check '     // &
                                 ' emission module documentation.'                                 
                CALL message( 'chem_emissions_setup', 'CM0446', 2, 2, 0, 6, 0 )

          END SELECT

       ENDIF

!
!-- Conversion factor to convert  kg/m**2/s to ppm/s

       DO  i = nxl, nxr
          DO  j = nys, nyn

!
!-- Derive Temperature from Potential Temperature

             tmp_temp(nzb:nzt+1,j,i) = pt(nzb:nzt+1,j,i) *                     &
                                       ( hyp(nzb:nzt+1) / p_0 )**rd_d_cp
 
!            
!-- We need to pass to cssws <- (ppm/s) * dz 
!-- Input is Nmole/(m^2*s)
!-- To go to ppm*dz multiply the input by (m**2/N)*dz
!-- (m**2/N)*dz == V/N 
!-- V/N = RT/P 

             conv_to_ratio(nzb:nzt+1,j,i) =  rgas_univ *                       &  ! J K-1 mol-1
                                             tmp_temp(nzb:nzt+1,j,i) /         &  ! K
                                             hyp(nzb:nzt+1)                       ! Pa

! (ecc) for reference
!                   m**3/Nmole               (J/mol)*K^-1           K                      Pa         
!             conv_to_ratio(nzb:nzt+1,j,i) = ( (Rgas * tmp_temp(nzb:nzt+1,j,i)) / ((hyp(nzb:nzt+1))) )  

          ENDDO
       ENDDO


! (ecc) moved initialization immediately after allocation
!
!-- Initialize

!       emis_distribution(:,nys:nyn,nxl:nxr,:) = 0.0_wp

  
!
!-- LOD 2 (PRE-PROCESSED MODE)

       IF  ( emiss_lod == 2 )  THEN

! for reference (ecc)
!       IF  ( TRIM( mode_emis ) == "PRE-PROCESSED" )  THEN

!
!-- Update time indices

          CALL get_date_time( 0.0_wp, second_of_day=time_utc_init )
          CALL get_date_time( MAX( 0.0_wp, time_since_reference_point ),       &
                              hour=hour_of_day )

          days_since_reference_point = INT( FLOOR( (                            &
                    time_utc_init + MAX( 0.0_wp, time_since_reference_point ) ) &
                                                   / seconds_per_day ) )
!*******************************************************************************
! MONA:
          index_hh = INT( FLOOR( MAX( 0.0_wp, time_since_reference_point ) / 3600.0 ) )
          write(9,*) 'index_hh = ', index_hh, time_since_reference_point
          write(9,*) dt_emission <= time_since_reference_point
          write(9,*) findloc( dt_emission <= time_since_reference_point, .TRUE. )
          write(9,*) maxval( findloc( dt_emission <= time_since_reference_point, .TRUE. ) )
          flush(9)
!          index_hh = days_since_reference_point * hours_per_day + hour_of_day
!*******************************************************************************

!
!-- LOD 1 (DEFAULT MODE)

        ELSEIF  ( emiss_lod == 1 )  THEN

! for reference (ecc)
!       ELSEIF  ( TRIM( mode_emis ) == "DEFAULT" )  THEN

!
!-- Allocate array where to store temporary emission values

          IF  (  .NOT.  ALLOCATED(emis) ) ALLOCATE( emis(nys:nyn,nxl:nxr) )

!
!-- Allocate time factor per category

          ALLOCATE( time_factor(emt_att%ncat) )

!
!-- Read-in hourly emission time factor

          IF  ( TRIM(time_fac_type) == "HOUR" )  THEN

!
!-- Update time indices

             CALL get_date_time( MAX( time_since_reference_point, 0.0_wp ),    &
                                 day_of_year=day_of_year, hour=hour_of_day )
             index_hh = ( day_of_year - 1_iwp ) * hour_of_day

!
!-- Check if the index is less or equal to the temporal dimension of HOURLY emission files

             IF  ( index_hh <= SIZE( emt_att%hourly_emis_time_factor(1,:) ) )  THEN

!
!-- Read-in the correspondant time factor

                time_factor(:) = emt_att%hourly_emis_time_factor(:,index_hh+1)      

!
!-- Error check (time out of range)

             ELSE

                message_string = 'The "HOUR" time-factors in the DEFAULT mode '           // &
                              ' are not provided for each hour of the total simulation time'      
                CALL message( 'chem_emissions_setup', 'CM0448', 2, 2, 0, 6, 0 ) 

             ENDIF

!
!-- Read-in MDH emissions time factors

          ELSEIF  ( TRIM( time_fac_type ) == "MDH" )  THEN

!
!-- Update time indices
             CALL get_date_time( MAX( time_since_reference_point, 0.0_wp ), &
                                 month=month_of_year,                       &
                                 day=day_of_month,                          &
                                 hour=hour_of_day,                          &
                                 day_of_week=day_of_week     )
             index_mm = month_of_year
             index_dd = months_per_year + day_of_week
             SELECT CASE(TRIM(daytype_mdh))

                CASE ("workday")
                   index_hh = months_per_year + days_per_week + hour_of_day

                CASE ("weekend")
                   index_hh = months_per_year + days_per_week + hours_per_day + hour_of_day

                CASE ("holiday")
                   index_hh = months_per_year + days_per_week + 2*hours_per_day + hour_of_day

             END SELECT
!
!-- Check if the index is less or equal to the temporal dimension of MDH emission files

             IF  ( ( index_hh + index_dd + index_mm) <= SIZE( emt_att%mdh_emis_time_factor(1,:) ) )  THEN

!
!-- Read in corresponding time factor

                time_factor(:) = emt_att%mdh_emis_time_factor(:,index_mm) *    &
                                 emt_att%mdh_emis_time_factor(:,index_dd) *    &
                                 emt_att%mdh_emis_time_factor(:,index_hh+1)

!
!-- Error check (MDH time factor not provided)
     
             ELSE

                message_string = 'The "MDH" time-factors in the DEFAULT mode '           //  &
                              ' are not provided for each hour/day/month  of the total simulation time'      
                CALL message( 'chem_emissions_setup', 'CM0449', 2, 2, 0, 6, 0 )

             ENDIF  

!
!-- Error check (no time factor defined)

          ELSE
                 
             message_string = 'In the DEFAULT mode the time factor'           //  &
                              ' has to be defined in the NAMELIST'      
             CALL message( 'chem_emissions_setup', 'CM0450', 2, 2, 0, 6, 0 )
         
          ENDIF

!
!-- PARAMETERIZED MODE

       ELSEIF ( emiss_lod == 0 )  THEN


! for reference (ecc)
!       ELSEIF  ( TRIM( mode_emis ) == "PARAMETERIZED" )  THEN
        
!
!-- assign constant values of time factors, diurnal profile for traffic sector

          par_emis_time_factor( : ) =  (/ 0.009, 0.004, 0.004, 0.009, 0.029, 0.039,   &
                                          0.056, 0.053, 0.051, 0.051, 0.052, 0.055,   &
                                          0.059, 0.061, 0.064, 0.067, 0.069, 0.069,   &
                                          0.049, 0.039, 0.039, 0.029, 0.024, 0.019 /)
          
          IF  (  .NOT.  ALLOCATED (time_factor) )  ALLOCATE (time_factor(1))

!
!--       Get time-factor for specific hour
          CALL get_date_time( MAX( time_since_reference_point, 0.0_wp ),              &
                              hour=hour_of_day )

          index_hh = hour_of_day
          time_factor(1) = par_emis_time_factor(index_hh+1)

       ENDIF  ! emiss_lod


!
!--  Emission distribution calculation

!
!-- LOD 0 (PARAMETERIZED mode)

       IF  ( emiss_lod == 0 )  THEN

! for reference (ecc)
!       IF  ( TRIM( mode_emis ) == "PARAMETERIZED" )  THEN

          DO  ispec = 1, n_matched_vars

!
!-- Units are micromoles/m**2*day (or kilograms/m**2*day for PMs)

             emis_distribution(1,nys:nyn,nxl:nxr,ispec) =            &
                       surface_csflux(match_spec_input(ispec)) *     &
                       time_factor(1) * hour_to_s

          ENDDO  

          
!
!-- LOD 1 (DEFAULT mode)


       ELSEIF  ( emiss_lod == 1 )  THEN

! for referene (ecc)
!       ELSEIF  ( TRIM( mode_emis ) == "DEFAULT" )  THEN

!
!-- Allocate array for the emission value corresponding to a specific category and time factor

          ALLOCATE (delta_emis(nys:nyn,nxl:nxr))

!
!-- Cycle over categories

          DO  icat = 1, emt_att%ncat
 
!
!-- Cycle over Species:  n_matched_vars represents the number of species
!-- in common between the emission input data and the chemistry mechanism used

             DO  ispec = 1, n_matched_vars

                emis(nys:nyn,nxl:nxr) =                    &
                       emt(match_spec_input(ispec))%       &
                           default_emission_data(icat,nys+1:nyn+1,nxl+1:nxr+1)

!
!-- NO 

                IF  ( TRIM(spc_names(match_spec_model(ispec))) == "NO" )  THEN
                
                   delta_emis(nys:nyn,nxl:nxr) = emis(nys:nyn,nxl:nxr) *       &         ! kg/m2/s
                                                 time_factor(icat) *           &
                                                 emt_att%nox_comp(icat,1) *    &
                                                 conversion_factor * hours_per_day

                   emis_distribution(1,nys:nyn,nxl:nxr,ispec) =                &
                               emis_distribution(1,nys:nyn,nxl:nxr,ispec) +    &
                               delta_emis(nys:nyn,nxl:nxr)
!
!-- NO2

                ELSEIF  ( TRIM(spc_names(match_spec_model(ispec))) == "NO2" )  THEN

                   delta_emis(nys:nyn,nxl:nxr) = emis(nys:nyn,nxl:nxr) *       &         ! kg/m2/s
                                                 time_factor(icat) *           &
                                                 emt_att%nox_comp(icat,2) *    &
                                                 conversion_factor * hours_per_day

                   emis_distribution(1,nys:nyn,nxl:nxr,ispec) =                &
                               emis_distribution(1,nys:nyn,nxl:nxr,ispec) +    &
                               delta_emis(nys:nyn,nxl:nxr)
 
!
!-- SO2
                ELSEIF  ( TRIM(spc_names(match_spec_model(ispec))) == "SO2" )  THEN
                  
                   delta_emis(nys:nyn,nxl:nxr) = emis(nys:nyn,nxl:nxr) *       &         ! kg/m2/s
                                                 time_factor(icat) *           &
                                                 emt_att%sox_comp(icat,1) *    &
                                                 conversion_factor * hours_per_day

                   emis_distribution(1,nys:nyn,nxl:nxr,ispec) =                &
                               emis_distribution(1,nys:nyn,nxl:nxr,ispec) +    &
                               delta_emis(nys:nyn,nxl:nxr)

!
!-- SO4

                ELSEIF  ( TRIM(spc_names(match_spec_model(ispec))) == "SO4" )  THEN

                  
                   delta_emis(nys:nyn,nxl:nxr) = emis(nys:nyn,nxl:nxr) *       &         ! kg/m2/s
                                                 time_factor(icat) *           &
                                                 emt_att%sox_comp(icat,2) *    &
                                                 conversion_factor * hours_per_day

                   emis_distribution(1,nys:nyn,nxl:nxr,ispec) =                &
                               emis_distribution(1,nys:nyn,nxl:nxr,ispec) +    &
                               delta_emis(nys:nyn,nxl:nxr)
 

!
!-- PM1

                ELSEIF  ( TRIM(spc_names(match_spec_model(ispec))) == "PM1" )  THEN

                   DO  i_pm_comp = 1, SIZE( emt_att%pm_comp(1,:,1) )   ! cycle through components

                      delta_emis(nys:nyn,nxl:nxr) = emis(nys:nyn,nxl:nxr) *              &     ! kg/m2/s
                                                    time_factor(icat) *                  &
                                                    emt_att%pm_comp(icat,i_pm_comp,1) *  &
                                                    conversion_factor * hours_per_day 

                      emis_distribution(1,nys:nyn,nxl:nxr,ispec) =                       &
                               emis_distribution(1,nys:nyn,nxl:nxr,ispec) +              &
                               delta_emis(nys:nyn,nxl:nxr)

                   ENDDO

!
!-- PM2.5

                ELSEIF  ( TRIM(spc_names(match_spec_model(ispec))) == "PM25" )  THEN

                   DO  i_pm_comp = 1, SIZE( emt_att%pm_comp(1,:,2) )   ! cycle through components

                      delta_emis(nys:nyn,nxl:nxr) = emis(nys:nyn,nxl:nxr) *              &     ! kg/m2/s
                                                    time_factor(icat) *                  &
                                                    emt_att%pm_comp(icat,i_pm_comp,2) *  &
                                                    conversion_factor * hours_per_day 

                      emis_distribution(1,nys:nyn,nxl:nxr,ispec) =                       &
                               emis_distribution(1,nys:nyn,nxl:nxr,ispec) +              &
                               delta_emis(nys:nyn,nxl:nxr)
 
                   ENDDO

!
!-- PM10

                ELSEIF  ( TRIM(spc_names(match_spec_model(ispec))) == "PM10" )  THEN

                   DO  i_pm_comp = 1, SIZE( emt_att%pm_comp(1,:,3) )   ! cycle through components

                      delta_emis(nys:nyn,nxl:nxr) = emis(nys:nyn,nxl:nxr) *              &     ! kg/m2/s
                                                    time_factor(icat)     *              &
                                                    emt_att%pm_comp(icat,i_pm_comp,3) *  &
                                                    conversion_factor * hours_per_day 

                      emis_distribution(1,nys:nyn,nxl:nxr,ispec) =                       &
                               emis_distribution(1,nys:nyn,nxl:nxr,ispec) +              &
                               delta_emis(nys:nyn,nxl:nxr) 

                   ENDDO

!
!-- VOCs

                ELSEIF  ( SIZE( match_spec_voc_input ) > 0 )  THEN

                   DO  ivoc = 1, SIZE( match_spec_voc_input )          ! cycle through components

                      IF  ( TRIM(spc_names(match_spec_model(ispec))) ==                &
                            TRIM(emt_att%voc_name(ivoc)) )  THEN    

                         delta_emis(nys:nyn,nxl:nxr) = emis(nys:nyn,nxl:nxr) *           &
                                                       time_factor(icat) *               &
                                                       emt_att%voc_comp(icat,match_spec_voc_input(ivoc)) *   &
                                                       conversion_factor * hours_per_day

                         emis_distribution(1,nys:nyn,nxl:nxr,ispec) =                    &
                               emis_distribution(1,nys:nyn,nxl:nxr,ispec) +              &
                               delta_emis(nys:nyn,nxl:nxr) 

                      ENDIF                       

                   ENDDO
                
!
!-- any other species

                ELSE

                   delta_emis(nys:nyn,nxl:nxr) = emis(nys:nyn,nxl:nxr) *                 &
                                                 time_factor(icat) *                     &
                                                 conversion_factor * hours_per_day
 
                   emis_distribution(1,nys:nyn,nxl:nxr,ispec) =                          &
                               emis_distribution(1,nys:nyn,nxl:nxr,ispec) +              &
                               delta_emis(nys:nyn,nxl:nxr) 

                ENDIF  ! TRIM spc_names
                
                emis = 0 
                
             ENDDO
             
             delta_emis = 0 
          
          ENDDO

!
!-- LOD 2 (PRE-PROCESSED mode)

       ELSEIF  ( emiss_lod == 2 )  THEN

! for reference (ecc)
!       ELSEIF  ( TRIM( mode_emis ) == "PRE-PROCESSED" )  THEN

!
!-- Cycle over species: n_matched_vars represents the number of species
!-- in common between the emission input data and the chemistry mechanism used

          DO  ispec = 1, n_matched_vars  
 
! (ecc)   
             emis_distribution(1,nys:nyn,nxl:nxr,ispec) =                                &
                       emt(match_spec_input(ispec))%                                     &
                           preproc_emission_data(index_hh+1,1,nys+1:nyn+1,nxl+1:nxr+1) * &
                       conversion_factor


!             emis_distribution(1,nys:nyn,nxl:nxr,ispec) =                                &
!                       emt(match_spec_input(ispec))%                                     &
!                           preproc_emission_data(index_hh,1,:,:) *   &
!                       conversion_factor
          ENDDO

       ENDIF  ! emiss_lod

       
!
!-- Cycle to transform x,y coordinates to the one of surface_mod and to assign emission values to cssws

!
!-- LOD 0 (PARAMETERIZED mode)
!-- Units of inputs are micromoles/m2/s

       IF  ( emiss_lod == 0 )  THEN
! for reference (ecc)
!       IF  ( TRIM( mode_emis ) == "PARAMETERIZED" )  THEN

          IF  (street_type_f%from_file)  THEN

!
!-- Streets are lsm surfaces, hence, no usm surface treatment required. 
!-- However, urban surface may be initialized via default initialization
!-- in surface_mod, e.g. at horizontal urban walls that are at k == 0
!-- (building is lower than the first grid point). Hence, in order to 
!-- have only emissions at streets, set the surfaces emissions to zero
!-- at urban walls.
 
             IF  ( surf_usm_h%ns >=1 )  surf_usm_h%cssws = 0.0_wp

!
!-- Treat land-surfaces.

             DO  m = 1, surf_lsm_h%ns

                i = surf_lsm_h%i(m)
                j = surf_lsm_h%j(m)
                k = surf_lsm_h%k(m)

!
!-- set everything to zero then reassign according to street type

                surf_lsm_h%cssws(:,m) = 0.0_wp

                IF  ( street_type_f%var(j,i) >= main_street_id    .AND.                  &
                      street_type_f%var(j,i) < max_street_id   )  THEN

!
!-- Cycle over matched species

                   DO  ispec = 1, n_matched_vars

!
!-- PMs are already in kilograms 

                      IF  ( TRIM(spc_names(match_spec_model(ispec))) == "PM1"     .OR.  &
                            TRIM(spc_names(match_spec_model(ispec))) == "PM25"    .OR.  &
                            TRIM(spc_names(match_spec_model(ispec))) == "PM10" )  THEN

!
!-- kg/(m^2*s) * kg/m^3
                         surf_lsm_h%cssws(match_spec_model(ispec),m) =                   &
                                   emiss_factor_main(match_spec_input(ispec)) *          &
                                   emis_distribution(1,j,i,ispec) *                      &     ! kg/(m^2*s)
                                   rho_air(k)                                                  ! kg/m^3
                         
!
!-- Other Species
!-- Inputs are micromoles

                      ELSE

!   
!-- ppm/s *m *kg/m^3                
                         surf_lsm_h%cssws(match_spec_model(ispec),m) =                   &
                                   emiss_factor_main(match_spec_input(ispec)) *          &
                                   emis_distribution(1,j,i,ispec) *                      &     ! micromoles/(m^2*s)
                                   conv_to_ratio(k,j,i) *                                &     ! m^3/Nmole     
                                   rho_air(k)                                                  ! kg/m^3

                      ENDIF

                   ENDDO  ! ispec


                ELSEIF  ( street_type_f%var(j,i) >= side_street_id   .AND.               &
                          street_type_f%var(j,i) < main_street_id )  THEN

!
!-- Cycle over matched species

                   DO  ispec = 1, n_matched_vars

!
!-- PMs are already in kilograms

                      IF  ( TRIM(spc_names(match_spec_model(ispec))) == "PM1"      .OR.  &
                            TRIM(spc_names(match_spec_model(ispec))) == "PM25"     .OR.  &
                            TRIM(spc_names(match_spec_model(ispec))) == "PM10" )  THEN

!
!-- kg/(m^2*s) * kg/m^3
                         surf_lsm_h%cssws(match_spec_model(ispec),m) =                   &
                                   emiss_factor_side(match_spec_input(ispec)) *          &
                                   emis_distribution(1,j,i,ispec) *                      &     ! kg/(m^2*s)
                                   rho_air(k)                                                  ! kg/m^3   
!
!-- Other species
!-- Inputs are micromoles

                      ELSE

!   
!-- ppm/s *m *kg/m^3

                         surf_lsm_h%cssws(match_spec_model(ispec),m) =                   &
                                   emiss_factor_side(match_spec_input(ispec)) *          &
                                   emis_distribution(1,j,i,ispec) *                      &  ! micromoles/(m^2*s)
                                   conv_to_ratio(k,j,i) *                                &  ! m^3/Nmole       
                                   rho_air(k)                                               ! kg/m^3   

                      ENDIF

                   ENDDO  ! ispec

!
!-- If no street type is defined, then assign zero emission to all the species

! (ecc) moved to front (for reference)
!                ELSE
!
!                   surf_lsm_h%cssws(:,m) = 0.0_wp

                ENDIF  ! street type

             ENDDO  ! m

          ENDIF  ! street_type_f%from_file


!
!-- LOD 1 (DEFAULT) and LOD 2 (PRE-PROCESSED)


       ELSE    
       

          DO  ispec = 1, n_matched_vars 
                   
!
!--          Default surfaces

             DO  m = 1, surf_def_h(0)%ns

                i = surf_def_h(0)%i(m)
                j = surf_def_h(0)%j(m)

                IF  ( emis_distribution(1,j,i,ispec) > 0.0_wp )  THEN

!
!-- PMs
                   IF  ( TRIM(spc_names(match_spec_model(ispec))) == "PM1"     .OR.    &
                         TRIM(spc_names(match_spec_model(ispec))) == "PM25"    .OR.    &
                         TRIM(spc_names(match_spec_model(ispec))) == "PM10" )  THEN
                
                      surf_def_h(0)%cssws(match_spec_model(ispec),m) =         &     ! kg/m2/s * kg/m3
                               emis_distribution(1,j,i,ispec)*                 &     ! kg/m2/s
                               rho_air(nzb)                                          ! kg/m^3  
 
                   ELSE

!
!-- VOCs
                      IF  ( len_index_voc > 0                                           .AND.      &
                            emt_att%species_name(match_spec_input(ispec)) == "VOC" )  THEN

                         surf_def_h(0)%cssws(match_spec_model(ispec),m) =      &  ! ppm/s * m * kg/m3
                               emis_distribution(1,j,i,ispec) *                &  ! mole/m2/s
                               conv_to_ratio(nzb,j,i) *                        &  ! m^3/mole
                               ratio2ppm *                                     &  ! ppm
                               rho_air(nzb)                                       ! kg/m^3  


!
!-- Other species

                      ELSE

                         surf_def_h(0)%cssws(match_spec_model(ispec),m) =      &   ! ppm/s * m * kg/m3
                               emis_distribution(1,j,i,ispec) *                &   ! kg/m2/s
                               ( 1.0_wp / emt_att%xm(ispec) ) *                &   ! mole/kg 
                               conv_to_ratio(nzb,j,i) *                        &   ! m^3/mole  
                               ratio2ppm *                                     &   ! ppm 
                               rho_air(nzb)                                        !  kg/m^3  
  
                      ENDIF  ! VOC

                   ENDIF  ! PM

                ENDIF  ! emis_distribution > 0

             ENDDO  ! m
         
!
!-- LSM surfaces

             DO  m = 1, surf_lsm_h%ns

                i = surf_lsm_h%i(m)
                j = surf_lsm_h%j(m)
                k = surf_lsm_h%k(m)

                IF  ( emis_distribution(1,j,i,ispec) > 0.0_wp )  THEN

!
!-- PMs
                   IF  ( TRIM(spc_names(match_spec_model(ispec))) == "PM1"     .OR.    &
                         TRIM(spc_names(match_spec_model(ispec))) == "PM25"    .OR.    &
                         TRIM(spc_names(match_spec_model(ispec))) == "PM10" )  THEN

                      surf_lsm_h%cssws(match_spec_model(ispec),m) =            &    ! kg/m2/s * kg/m3
                               emis_distribution(1,j,i,ispec) *                &    ! kg/m2/s
                               rho_air(k)                                           ! kg/m^3
 
                   ELSE

!
!-- VOCs

                      IF  ( len_index_voc > 0                                           .AND.      &
                            emt_att%species_name(match_spec_input(ispec)) == "VOC" )  THEN

                         surf_lsm_h%cssws(match_spec_model(ispec),m) =         &   ! ppm/s * m * kg/m3
                               emis_distribution(1,j,i,ispec) *                &   ! mole/m2/s  
                               conv_to_ratio(k,j,i) *                          &   ! m^3/mole
                               ratio2ppm *                                     &   ! ppm
                               rho_air(k)                                          ! kg/m^3  

!
!-- Other species

                      ELSE

                         surf_lsm_h%cssws(match_spec_model(ispec),m) =         &   ! ppm/s * m * kg/m3
                               emis_distribution(1,j,i,ispec) *                &   ! kg/m2/s
                               ( 1.0_wp / emt_att%xm(ispec) ) *                &   ! mole/kg
                               conv_to_ratio(k,j,i) *                          &   ! m^3/mole
                               ratio2ppm *                                     &   ! ppm
                               rho_air(k)                                          ! kg/m^3     
                                                   
                      ENDIF  ! VOC

                   ENDIF  ! PM

                ENDIF  ! emis_distribution

             ENDDO  ! m

!
!-- USM surfaces

             DO  m = 1, surf_usm_h%ns

                i = surf_usm_h%i(m)
                j = surf_usm_h%j(m)
                k = surf_usm_h%k(m)

                IF  ( emis_distribution(1,j,i,ispec) > 0.0_wp )  THEN

!
!-- PMs
                   IF  ( TRIM(spc_names(match_spec_model(ispec))) == "PM1"     .OR.    &
                         TRIM(spc_names(match_spec_model(ispec))) == "PM25"    .OR.    &
                         TRIM(spc_names(match_spec_model(ispec))) == "PM10" )  THEN
                   
                      surf_usm_h%cssws(match_spec_model(ispec),m) =            &    ! kg/m2/s * kg/m3
                               emis_distribution(1,j,i,ispec)*                 &    ! kg/m2/s
                               rho_air(k)                                           ! kg/m^3

                   ELSE
                      
!
!-- VOCs
                      IF  ( len_index_voc > 0                                         .AND.      &
                            emt_att%species_name(match_spec_input(ispec)) == "VOC" )  THEN

                         surf_usm_h%cssws(match_spec_model(ispec),m) =         &   ! ppm/s * m * kg/m3
                               emis_distribution(1,j,i,ispec) *                &   ! m2/s
                               conv_to_ratio(k,j,i) *                          &   ! m^3/mole
                               ratio2ppm *                                     &   ! ppm
                               rho_air(k)                                          ! kg/m^3   

!
!-- Other species
                      ELSE

                         surf_usm_h%cssws(match_spec_model(ispec),m) =         &   ! ppm/s * m * kg/m3
                               emis_distribution(1,j,i,ispec) *                &   ! kg/m2/s
                               ( 1.0_wp / emt_att%xm(ispec) ) *                &   ! mole/kg
                               conv_to_ratio(k,j,i) *                          &   ! m^3/mole
                               ratio2ppm*                                      &   ! ppm
                               rho_air(k)                                          ! kg/m^3   


                      ENDIF  ! VOC

                   ENDIF  ! PM

                ENDIF  ! emis_distribution
 
             ENDDO  ! m

          ENDDO

       ENDIF

!
!-- Deallocate time_factor in case of DEFAULT mode) 

       IF  ( ALLOCATED (time_factor) )  DEALLOCATE (time_factor)

   ENDIF

 END SUBROUTINE chem_emissions_setup

 END MODULE chem_emissions_mod
