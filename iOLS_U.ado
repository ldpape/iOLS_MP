** 28/06/2021 : Corrected error on Diagonal Matrix of Weights using "cross".
** 14/12/2021 : Corrected convergence criteria with "( . )"
** 14/12/2021 : Changed Convergence Criteria from Absolute change to Relative Change
** 14/12/2021 : Added a quietly after "preserve" 
** 14/12/2021 : Changed the constant calculation to avoid numerical log(0).
** 21/12/2021 : Updated to matrix form for speed and options to control convergence.
** 04/01/2021 : Add additional stopping criteria + return of the constant alpha.
** 05/01/2021 : Created iOLS_U

cap program drop iOLS_U
program define iOLS_U, eclass 
	syntax [anything] [if] [in] [aweight pweight fweight iweight] [, DELta(real 0) Robust LIMit(real 0.00001) MAXimum(real 1000) CLuster(varlist numeric)]
	marksample touse
	preserve
	quietly keep if `touse'
	if  "`robust'" !="" {
		local opt1  = "`robust' "
	}
	if "`cluster'" !="" {
		local opt2 = "vce(cluster `cluster') "
	}
	local option = "`opt1'`opt2'"
	local list_var `anything'
	* get depvar and indepvar
	gettoken depvar list_var : list_var
	gettoken indepvar list_var : list_var, p("(")
	tempvar y_tild 
	quietly gen `y_tild' = log(`depvar' + 1)
	tempvar cste
	gen `cste' = 1
	** drop collinear variables
    _rmcoll `indepvar' `cste', forcedrop 
	local var_list `endog' `r(varlist)' `cste'  
	*** Initialisation de la boucle
	mata : X=.
	mata : y_tilde =.
	mata : y =.
	mata : st_view(X,.,"`var_list'")
	mata : st_view(y_tilde,.,"`y_tild'")
	mata : st_view(y,.,"`depvar'")
	mata : invXX = invsym(cross(X,X))
	mata : beta_initial = invXX*cross(X,y_tilde)
	mata : beta_t_1 = beta_initial // needed to initialize
	mata : beta_t_2 = beta_initial // needed to initialize
	** initialize delta path
	mata: scale_scalar = max(y:*exp(-X*beta_initial))
	if `delta'==0{
	mata : delta= 0.01*scale_scalar // max((`delta', 0.5*scale_scalar ))
	}
	else{
	mata: delta = `delta'
	}
	mata : q_hat_m0 = 0
	local k = 1
	local eps = 1000	
	mata: q_hat = J(`maximum', 1, .)
	*** ItÃ©rations iOLS
	_dots 0
	while ( (`k' < `maximum') & (`eps' > `limit') ) {
	*mata: alpha = log(mean(y:*exp(-X[.,1..(cols(X)-1)]*beta_initial[1..(cols(X)-1),1])))
	*mata : beta_initial[(cols(X)),1] = alpha
	mata: xb_hat = X*beta_initial
		* Update d'un nouveau y_tild et regression avec le nouvel y_tild
	mata: y_tilde = log(y + delta*exp(xb_hat)) :- (log(delta :+ y:*exp(-xb_hat)) :- ((y:*exp(-xb_hat) :- 1):/(1:+delta)))
		* 2SLS 
	mata: beta_new = invXX*cross(X,y_tilde)
		* Difference entre les anciens betas et les nouveaux betas
	mata: criteria = mean(abs(beta_initial - beta_new):^(2))
mata: st_numscalar("eps", criteria)
mata: st_local("eps", strofreal(criteria))
		* safeguard for convergence.
	if `k'==`maximum'{
		  di "There has been no convergence so far: increase the number of iterations."  
	}
	if `k'>4{
	mata: q_hat[`k',1] = mean(log( abs(beta_new-beta_initial):/abs(beta_initial-beta_t_2)):/log(abs(beta_initial-beta_t_2):/abs(beta_t_2-beta_t_3)))	
	mata: check_3 = abs(mean(q_hat)-1)
	if mod(`k'-4,50)==0{
*	mata: if (hasmissing(q_hat[((`k'-49)..`k'),.])==0) q_hat_m=mm_median(q_hat[((`k'-49)..`k'),.] ,1) ;;
	mata: q_hat_m=mm_median(select(q_hat[((`k'-49)..`k'),.],rowmissing(q_hat[((`k'-49)..`k'),.]):==0) ,1)
	mata: check_1 = abs(q_hat_m - q_hat_m0)
	mata: check_2 = abs(q_hat_m-1)
	mata: q_hat_m0 = q_hat_m ;;
	mata: st_numscalar("check_1", check_1)
	mata: st_local("check_1", strofreal(check_1))
	mata: st_numscalar("check_2", check_2)
	mata: st_local("check_2", strofreal(check_2))
	mata: st_numscalar("check_3", check_3)
	mata: st_local("check_3", strofreal(check_3))
		if ((`check_1'<1e-4)&(`check_2'>1e-2)) {
di "delta is too small to achieve convergence -- updating to larger value"
	*local k = `maximum'
	mata : delta = delta + 0.66*scale_scalar
		}
		if ((`check_3'>0.5) & (`k'>50)) {
*	local k = `maximum'
di "q_hat too far from 1: updating to larger value"
mata : delta = delta+ 0.66*scale_scalar
		}
					  }
	}
	if `k'>2 { // keep in memory the previous beta_hat for q_hat 
	mata:   beta_t_3 = beta_t_2
	mata:   beta_t_2 = beta_initial
	}
mata: beta_initial = beta_new
	local k = `k'+1
	_dots `k' 0
	}

	*** Calcul de la bonne matrice de variance-covariance
	* Calcul du "bon" rÃ©sidu
	mata: xb_hat = X*beta_new
	mata: y_tilde = log(y + delta*exp(xb_hat)) :- (log(delta :+ y:*exp(-xb_hat)) :- ((y:*exp(-xb_hat) :- 1):/(1:+delta)))
	mata: ui = y:*exp(-xb_hat)
	mata: ui = ui:/(delta :+ 1)
	* Retour en Stata 
	cap drop y_tild 
	quietly mata: st_addvar("double", "y_tild")
	mata: st_store(.,"y_tild",y_tilde)
	quietly reg y_tild `r(varlist)' [`weight'`exp'] if `touse', `option'
	matrix beta_final = e(b) // 	mata: st_matrix("beta_final", beta_new)
	matrix Sigma = e(V)
	mata : Sigma_hat = st_matrix("Sigma")
	mata : Sigma_0 = cross(X,X)*Sigma_hat*cross(X,X)
	mata : invXpIWX = invsym(cross(X, ui, X))  
	mata : Sigma_tild = invXpIWX*Sigma_0*invXpIWX
    mata: st_matrix("Sigma_tild", Sigma_tild)
	mata: st_numscalar("delta", delta)
mata: st_local("delta", strofreal(delta))
	*** Stocker les resultats dans une matrice
	local names : colnames e(b)
	local nbvar : word count `names'
	mat rownames Sigma_tild = `names' 
    mat colnames Sigma_tild = `names' 
    ereturn post beta_final Sigma_tild , obs(`=e(N)') depname(`depvar') esample(`touse')  dof(`=r(df r)') 
	restore 
ereturn scalar delta = `delta'
ereturn  scalar eps =   `eps'
ereturn  scalar niter =  `k'
ereturn local cmd "iOLS_U"
ereturn local vcetype `option'
di in gr _col(55) "Number of obs = " in ye %8.0f e(N)
di in gr _col(55) "delta of conv = " in ye %8.0f delta
ereturn display

end

