    ! S8 (sigma8 * Omegam**pwr) measurement constraint

    module S8
    use CosmologyTypes
    use CosmoTheory
    use Likelihood_Cosmology
    implicit none
    private

    type, extends(TCosmoCalcLikelihood) :: S8Likelihood
        real(mcp) :: sm8, sm8_err, sm8_omt, sm8_pwr 
    contains
    procedure :: LogLike => S8_LnLike
    end type S8Likelihood


    public S8Likelihood, S8Likelihood_Add
    contains

    subroutine S8Likelihood_Add(LikeList, Ini)
    class(TLikelihoodList) :: LikeList
    class(TSettingIni) :: ini
    Type(S8Likelihood), pointer :: this

    if (Ini%Read_Logical('use_S8',.false.)) then
        allocate(this)
        this%LikelihoodType = 'S8'
        this%name= Ini%Read_String('S8_name')
        this%sm8 = Ini%Read_Double('S8_sm8')
        this%sm8_err = Ini%Read_Double('S8_sm8_err')
        this%sm8_omt = Ini%Read_Double('S8_sm8_omt')
        this%sm8_pwr = Ini%Read_Double('S8_sm8_pwr')
        call LikeList%Add(this)
    end if

    end subroutine S8Likelihood_Add

    real(mcp) function S8_LnLike(this, CMB, Theory, DataParams)
    Class(S8Likelihood) :: this
    Class(CMBParams) CMB
    Class(TCosmoTheoryPredictions), target :: Theory
    real(mcp) :: DataParams(:)
    real(mcp) :: theoryval

    if(Theory%sigma_8==0) &
                call MpiStop('ERROR: sigma_8 have not been computed. S8_LnLike')

    theoryval = Theory%sigma_8*( (CMB%omb+CMB%omdm)/this%sm8_omt)**this%sm8_pwr
   
    S8_LnLike = (theoryval - this%sm8)**2/(2*this%sm8_err**2)

    end function  S8_LnLike

    end module S8
