!----------------------------------------------------------------------------------------------------------------
!N
!  Action:   QM Harmonic Oscillator  N_t*eps=beta   here    eps = beta/N_t
!    
!----------------------------------------------------------------------------------------------------------------     
	PROGRAM MAIN 
        USE MUMBERS      
                                         
	IMPLICIT NONE
	   include 'mpif.h'      
!----------------------------Static/Helper_Variables-----------------------------------------------------------------------
        integer, parameter             :: d=1          		                                  ! one time dimension
        integer                        :: N_t                                                     ! Linear system size in time direction
	double precision               :: eps,beta,osc_freq,osc_mass,a_0                 		 	    
        double precision               :: integral,uplim,lowlim                       
        double precision               :: RF2(0:100000)                 
        double precision               :: step,scan_step,s_step,p_step,t_step,i_t,i_p,i_s               
        integer                        :: charge,j_1,prnt_number         
        character*256                  :: results_name,cnf_name,stats_name,fname,error_file,fef       
!------------------------Configuration_Variables-----------------------------------------------------------------------
        integer, allocatable           :: site_charges(:),bonds_t(:)                 ! arrays        
        integer                        :: t_m1,t_m2,shifter,wiggler             ! coordinates ...
!------------------------Statistical_Variables--------------------------------------------------------------------------
        ! double precision, allocatable  :: 
         double precision               :: Z,xx,Q_xx
         double precision               :: amax, tmax, amin, tmin
         integer                        :: stats_mode,config_mode,maxcharge
         integer, parameter             :: prt_m=100000                                          
         double precision               :: Q_1_prntouts(prt_m)
           
!------------------INITIALIZE MPI------------------------------------------------------------
        integer :: rank, comm_size, ierr, outproc, status(MPI_STATUS_SIZE)

        call MPI_INIT( ierr )
        call MPI_COMM_RANK( MPI_COMM_WORLD, rank, ierr )
        call MPI_COMM_SIZE( MPI_COMM_WORLD, comm_size, ierr )        
!------------INPUTS from PAR file-----------------------------------------------------------------
         CALL READ_INPUT_FILE         
!-------------------Filenames for MPI----------------------------------------------------------         
         if (N_t>99) then
           if (rank>9) then
             write(results_name,'(A,i3,A,i2,A,F3.1,A)') "output2/",N_t,"_**rank",rank,"_beta",beta,".dat" 
           else
             write(results_name,'(A,i3,A,i1,A,F3.1,A)') "output2/",N_t,"_**rank",rank,"_beta",beta,".dat"
           end if
         else if (N_t > 9) then
            if (rank>9) then
         	write(results_name,'(A,i2,A,i2,A,F3.1,A)') "output2/",N_t,"_**rank",rank,"_beta",beta,".dat" 
            else
            	write(results_name,'(A,i2,A,i1,A,F3.1,A)') "output2/",N_t,"_**rank",rank,"_beta",beta,".dat"
            end if
         else
          if (rank>9) then
             write(results_name,'(A,i1,A,i2,A,F3.1,A)') "output2/",N_t,"_**rank",rank,"_beta",beta,".dat"  
          else
             write(results_name,'(A,i1,A,i1,A,F3.1,A)') "output2/",N_t,"_**rank",rank,"_beta",beta,".dat"
          end if
         end if                                                 
 
        if(N_t>99) then                               !  formatting this will break if numprocs > 99
            if(rank>9) then                 
                write(cnf_name,1343) rank,'_N_t', N_t
                write(stats_name,3343) rank,'_N_t', N_t
               write(fname,4343) rank,'_N_t', N_t
            else                
                write(cnf_name,1344) rank,'_N_t', N_t
                write(stats_name,3344) rank,'_N_t', N_t
                write(fname,4344) rank,'_N_t', N_t
            end if
        else
            if(rank>9) then                 
               write(cnf_name,1346) rank,'_N_t', N_t
                write(stats_name,3346) rank,'_N_t', N_t
                write(fname,4346) rank,'_N_t', N_t
            else                
                write(cnf_name,1347) rank,'_N_t', N_t
                write(stats_name,3347) rank,'_N_t', N_t
                write(fname,4347) rank,'_N_t', N_t
            end if                
        endif
        
        1343    format('config/cnf_',i2,A,i3)
        3343    format('config/stats_',i2,A,i3)
        4343    format('config/rand_state_',i2,A,i3)                 
        1344    format('config/cnf_',i1,A,i3)
        3344    format('config/stats_',i1,A,i3)
        4344    format('config/rand_state_',i1,A,i3)                
        1346    format('config/cnf_',i2,A,i2)
        3346    format('config/stats_',i2,A,i2)
        4346    format('config/rand_state_',i2,A,i2)        
        1347    format('config/cnf_',i1,A,i2)
        3347    format('config/stats_',i1,A,i2)
        4347    format('config/rand_state_',i1,A,i2)               
!--------------Scan Parameter (for different MPI_processes)--------------------------------------------
        osc_freq = osc_freq + scan_step*rank                        
!------------calculate site factor integrals-------------------------------------------------------------------------------            
	      
        DO  j_1 = 0,100000   
             RF2(j_1) = (2*j_1+1)/(2*a_0)
        END DO
                                                                                                    
      if (config_mode == 1) then      
         CALL READ_CONFIG         
      else      
         CALL INIT_CNF         
      end if
      
      if (stats_mode == 1) then      
         CALL READ_STATS         
      else          
         CALL INIT_STAT    
      end if
       
      if (config_mode == 0 .AND. stats_mode == 0) then
	 DO  	                                ! annealling
	     if ( rndm() < 0.5d0 ) then		    
			shifter    = 0                                 
			wiggler = t_m1						
	     else		    
			shifter = 1                                  
			wiggler = t_m2
	     end if	    	         	        
	            CALL SHIFT         	  				
		        i_t = i_t + 1.d0		               			                   
			    if (i_t >= t_step) then		    
			        exit			        
			    end if	                                                    		
          END DO		     			    		    				                       		    	             
       end if    
!------------------------------main loop---------------------------------------------------------------------                                                                                                                                              
	    DO  
	    
	     if ( rndm() < 0.5d0 ) then		    
			shifter = 0                                 
			wiggler = t_m1
	     else		    
			shifter = 1                                  
			wiggler = t_m2
	     end if
	    	                 
	            CALL SHIFT         	  
                    CALL MEASURE
	            step = step + 1.d0	          
	            i_p = i_p + 1.d0
	            i_s = i_s + 1.d0	               
                                   
              if (i_p >= p_step) then                            !=> write the results, including errors                   
                  i_p = 0.d0                                    
                  CALL PRNT_RESULTS                 
              end if
              
              if (i_s  >= s_step)  then                          !=> save the configuration and statistics            
                  i_s = 0.d0                
                  CALL SAVE_CONFIG 
                  CALL SAVE_STATS                     
              end if  
                         		
        END DO
               
        call MPI_FINALIZE(ierr)
                       
    CONTAINS
!---------------------------------------------------------------------------------------------

!------------------------NUMERICAL-INTEGRATION------------------------------------------------
  ! not needed
!---------------------------------------------------------------------------------------------	
	    SUBROUTINE READ_INPUT_FILE 
	    
        OPEN(1, FILE='par')
        
        READ(1,*) N_t        
	READ(1,*) osc_mass
	READ(1,*) osc_freq
        READ(1,*) beta           
        READ(1,*) scan_step
        READ(1,*) t_step                       ! thermolization
        READ(1,*) p_step                       ! print
        READ(1,*) s_step                       ! save
        READ(1,*) config_mode
        READ(1,*) stats_mode
        
        ALLOCATE(bonds_t(0:N_t-1))        
        ALLOCATE(site_charges(0:N_t-1))
        
        eps = beta/N_t
        a_0 = (0.5d0*eps*eps*osc_freq*osc_freq + 1)
        
        t_step = t_step*d*dble(N_t)
        p_step = p_step*d*dble(N_t)
        s_step = s_step*d*dble(N_t)
            
        END SUBROUTINE READ_INPUT_FILE

!---------------------------------------------------------------------------------
! algorithm: randomly choose to shift either ira or masha
 
	     SUBROUTINE SHIFT
         
        integer :: direction,draw_or_erase,target_t,bond_number        		                                                                
        double precision :: acc_r
        
         direction = RN(2*d)                               ! select a direction at random
         
         if (direction == 2) then
             direction = -1
         end if
         
         draw_or_erase = RN(2)
         
         if (draw_or_erase == 2) then
             draw_or_erase = -1
         end if
          
        if (direction == 1) then                         ! PBC and bond number lookup
            target_t = wiggler + 1
            
            if (target_t == N_t) then
                target_t = 0
            end if			
            bond_number = bonds_t(wiggler)
                             
        else if (direction == -1) then
            target_t = wiggler - 1

            if (target_t == -1) then
                target_t = N_t-1
            end if
            bond_number = bonds_t(target_t)
        
        end if
        
!---------ACCEPTANCE RATIOS--------------------------------------------------------------------------------
        
        if (draw_or_erase == 1) then        
                acc_r = 1.d0/(1.d0 + bond_number)*RF2(site_charges(target_t))       
        else if (draw_or_erase == -1) then        
                acc_r = (bond_number*1.d0)*1.d0/RF2(site_charges(wiggler)-1)                     
        end if
                
!------------DECISION-TO-UPDATE/UPDATE------------------------------------------------------------------

        if (rndm() < acc_r) then
        	
        	!xsq = xsq - RF2(site_charges(target_t))
        	!xsq = xsq - RF2(site_charges(wiggler))
                                           
            if (direction == 1) then
                bonds_t(wiggler) = bonds_t(wiggler) + draw_or_erase          
            else if (direction == -1) then
                bonds_t(target_t)   = bonds_t(target_t) + draw_or_erase
           end if 
            
            if (draw_or_erase == 1) then               
                site_charges(target_t) = site_charges(target_t) + 1                            
            else if (draw_or_erase == -1) then            
                site_charges(wiggler) = site_charges(wiggler) - 1                          
            end if
            
            !xsq = xsq + RF2(site_charges(target_t))     
            !xsq = xsq + RF2(site_charges(wiggler))
            
            if (shifter == 0) then
	        t_m1 = target_t
                
            else if (shifter == 1) then
	        t_m2 = target_t		
	    end if
		    		    
          end if
        
        END SUBROUTINE SHIFT
!---------------------------------------------------------------------------------------------------------------------
        SUBROUTINE INIT_STAT       
        ! initialize statistical quantities
        
        Z   = 0.d0 
        xx = 0.d0
        Q_xx = 0.d0
        
        maxcharge = maxval(site_charges)

        step = 0.d0
        i_p  = 0.d0
        i_s  = 0.d0
        i_t  = 0.d0
        
        prnt_number = 0

        Q_1_prntouts = 0.d0

              
        END SUBROUTINE INIT_STAT
!-------------------------------------------------------------------------
!		Initializing a new configuration
!-------------------------------------------------------------------------
	    SUBROUTINE INIT_CNF
	    
        integer :: ij,kl,jj
        
!       Initializes configuration, sets all the configuration variables to 0  
                                          
         bonds_t = 0
 
         site_charges = 0
                       
         t_m1 = RN(N_t)                                               
         t_m1 = t_m1 -1                  
         t_m2 = t_m1
         site_charges(t_m1) =  1
                                                 
!       Initializing random number generator  
                                                                           
         ij = 1802 + 18*rank
         kl = 9373 + 17*rank
         CALL SRANMAR(ij,kl)     
        
        END SUBROUTINE INIT_CNF
!-------------------------------------------------------------------------------------------------------------
        SUBROUTINE MEASURE
                                   
           if ( t_m1 == t_m2 ) then                          	   
           	   Z = Z + 1.d0/RF2(site_charges(t_m1)-1)                                             !estimator for the partition function    
           	   xx = xx + 1.d0          	            	    
           end if
                                                                   
	    END SUBROUTINE MEASURE	    
!-------------------OUTPUT------------------------------------------
        SUBROUTINE PRNT_RESULTS
                
        double precision :: error_set_1(0:6)
                                                          		
           prnt_number = prnt_number + 1
		
	   if (maxval(site_charges) > maxcharge) then
           	   maxcharge = maxval(site_charges)
           end if         
           
        Q_1_prntouts(prnt_number) = xx/Z                                    
        call STAT(Q_1_prntouts,prnt_number,amax,tmax,amin,tmin)                          
        error_set_1(0) = prnt_number
        error_set_1(1) = Q_1_prntouts(prnt_number)
        error_set_1(2) = amax-amin                                        
        error_set_1(3) = amax                                             
        error_set_1(4) = tmax
        error_set_1(5) = amin
        error_set_1(6) = tmin

       
               open(53+rank,file = results_name)
         
             write(53+rank,*) ' '
             write(53+rank,'(A,i2,A,i3,A,F7.3,A,F7.3,A,ES10.3E2)')'rank = ',rank, ' size N_t =',N_t ,' beta = ', beta,'  m = ',osc_mass,'  mc run steps', step
             write(53+rank,*) ' '
             write(53+rank,*) ' w          <x^2>   '
             write(53+rank,*) ' '             
             write(53+rank,'(F8.4,A,F10.3,A,F10.3)') osc_freq,"   ", eps*xx/Z 
             write(53+rank,*) ' '  
             write(53+rank,'(A,i7,A,i4)') 'number of printouts: ', prnt_number, '   maximum charge = ', maxcharge
             write(53+rank,*) ' '
             write(53+rank,'(A,F12.5,A,F12.5)') ' <x^2> = ', error_set_1(1), '  max-min = ', error_set_1(2)
             write(53+rank,'(F12.5,A,F12.5)') error_set_1(3),' at -> ',error_set_1(4)
             write(53+rank,'(F12.5,A,F12.5)') error_set_1(5),' at -> ',error_set_1(6)
             write(53+rank,*) ' '

             	close(53+rank) 
          
	    END SUBROUTINE PRNT_RESULTS
!-----------------stat-------------------------------------------------------------------------------------------
      SUBROUTINE STAT(prnt_storage,n,amax,tmax,amin,tmin) 
!     Analyzing 3/4 print-out
	      
      double precision :: prnt_storage, amax, tmax, amin, tmin, aa
      integer n, i
      dimension prnt_storage(n)

      amax=-1.d200; amin=1.d200

      DO i=n/4+1, n;  aa=prnt_storage(i)
         if (aa > amax) then
            amax=aa; tmax=i
         end if
         if (aa < amin) then
            amin=aa; tmin=i
         end if
      END DO

      tmax=tmax/n; tmin=tmin/n
      END SUBROUTINE STAT	    
!------------SAVE CONFIGURATION--------------------------------------------------------------------
       
       SUBROUTINE SAVE_CONFIG
                  
       open(14, file=cnf_name)  
                               
       write(14,*) t_m1   
       write(14,*) t_m2
     
       write(14,*) bonds_t       

  
       write(14,*) site_charges
       
       close(14)                                                  
                                                     	                   						   
	   call save_rndm(fname)                     
	                                              ! (this uses file reference number 23)      
       END SUBROUTINE SAVE_CONFIG       
!--------------------Save Stats -----------------------------------------------------------------------------  
     SUBROUTINE SAVE_STATS
     
     open(15, file=stats_name)
     
     write(15,*) step  
     write(15,*) i_s
     write(15,*) i_p 

     write(15,*) Z

     write(15,*) Q_1_prntouts
     write(15,*) prnt_number
     write(15,*) maxcharge
     
     close(15)
         
     END SUBROUTINE SAVE_STATS
!----------------------------------------------------------------------------------------------------------------------
       SUBROUTINE READ_CONFIG
       
       open(14, file=cnf_name)  
                               
       read(14,*) t_m1   
       read(14,*) t_m2
     
       read(14,*) bonds_t       


       read(14,*) site_charges
       
       close(14)
       
       call read_rndm(fname)
                    
       END SUBROUTINE READ_CONFIG
!----------------------------------------------------------------------------------------------------------------------
       SUBROUTINE READ_STATS
       
     open(15, file=stats_name)
     
     read(15,*) step  
     read(15,*) i_s
     read(15,*) i_p  
     
     read(15,*) Z


     read(15,*) Q_1_prntouts
     read(15,*) prnt_number
     read(15,*) maxcharge
     
     close(15)
              
       END SUBROUTINE READ_STATS
!----------------------------------------------------------------------------------------------------------------------                   	
      END PROGRAM MAIN
      