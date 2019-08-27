/*! \file
 *  \brief Plaquette gauge action as sum of characters
 */

#include "chromabase.h"
#include "actions/gauge/gaugeacts/rmc_coupled_gaugeact.h"
#include "actions/gauge/gaugeacts/gaugeact_factory.h"
#include "actions/gauge/gaugestates/gauge_createstate_factory.h"
#include "actions/gauge/gaugestates/gauge_createstate_aggregate.h"
#include "io/aniso_io.h"

namespace Chroma
{
 
  namespace RMCCoupledGaugeActEnv 
  { 
    namespace
    {
      GaugeAction< multi1d<LatticeColorMatrix>, 
		   multi1d<LatticeColorMatrix> >* createGaugeAct(XMLReader& xml, 
								 const std::string& path) 
      {
	return new RMCCoupledGaugeAct(CreateGaugeStateEnv::reader(xml, path),
			    Params(xml, path));
      }
      
      const std::string name = "RMC_COUPLED_GAUGEACT";

      //! Local registration flag
      static bool registered = false;
    }

    //! Register all the factories
    bool registerAll() 
    {
      bool success = true; 
      if (! registered)
      {
	success &= TheGaugeActFactory::Instance().registerObject(name, createGaugeAct);
	registered = true;
      }
      return success;
    }


    Params::Params(XMLReader& xml_in, const std::string& path) 
    {
      XMLReader paramtop(xml_in, path);

      try 
      {	
	read(paramtop, "beta", beta);
	read(paramtop, "alpha", alpha);
	read(paramtop, "gamma", gamma);

	RMCGamma = zero;
	shadow_gauge_field = zero;
	read(paramtop, "shadow_gauge_id",shadow_gauge_id);
	
	std::string temp;
	read(paramtop, "file_format", temp);

	if(temp.compare("SERIAL") == 0) {
		file_format = QDPIO_SERIAL;
	}
	else if(temp.compare("PARALLEL") == 0) {
		file_format = QDPIO_PARALLEL;
	}
	else {
		QDPIO::cerr << "Error reading file format; please put either SERIAL or PARALLEL: " << std::endl;
		QDP_abort(1);
	}

	curr_num = 1;
	read(paramtop, "freeze_steps", freeze_steps);
	read(paramtop, "freeze_check", freeze_check);

      }
      catch( const std::string& e ) { 
	QDPIO::cerr << "Error reading XML: " <<  e << std::endl;
	QDP_abort(1);
      }

    }


    //! Compute the action
    Double RMCCoupledGaugeAct::S(const Handle< GaugeState<P,Q> >& state) const
    {
      // Action at the site level
      multi2d<LatticeColorMatrix> plq;
      this->siteAction(plq, state);

      // Total action
      Double act_F = zero;
      
      Double three = Nc;       

      for(int mu=1; mu < Nd; ++mu)
      {
	for(int nu=0; nu < mu; ++nu)
	{
	  LatticeComplex plaq=trace(plq[mu][nu]);
	  // Sum over plaquettes
	  act_F += sum(real(plaq)-three);		 
	}
      }

      // Normalize
      Real act = -(param.beta / Real(Nc)) * act_F;

      return act;
    }
 

    //! Compute the plaquette
    void RMCCoupledGaugeAct::siteAction(multi2d<LatticeColorMatrix>& plq, const Handle< GaugeState<P,Q> >& state) const
    {
      START_CODE();

      // Initialize
      plq.resize(Nd,Nd);
      plq = zero;

      // Handle< const GaugeState<P,Q> > u_bc(createState(u));
      // Apply boundaries
      const multi1d<LatticeColorMatrix>& u = state->getLinks();

      // Compute the average plaquettes
      for(int mu=1; mu < Nd; ++mu)
      {
	for(int nu=0; nu < mu; ++nu)
	{
	  plq[mu][nu] += (u[mu]*shift(u[nu],FORWARD,mu)*adj(shift(u[mu],FORWARD,nu))*adj(u[nu])) ; 

	  // Keep a copy
	  plq[nu][mu] = plq[mu][nu];
	}
      }

      END_CODE();
    }
 

    //! Compute staple
    /*! Default version. Derived class should override this if needed. */
    void RMCCoupledGaugeAct::staple(LatticeColorMatrix& result,
			  const Handle< GaugeState<P,Q> >& state,
			  int mu, int cb) const
    {
     START_CODE();

    const multi1d<LatticeColorMatrix>& u = state->getLinks();
    
    result = zero;
    LatticeColorMatrix tmp1, tmp2;
    LatticeColorMatrix u_nu_mu;

    for(int nu=0; nu < Nd; ++nu) 
    {
      if( nu == mu ) continue;
      
      u_nu_mu = shift(u[nu],FORWARD,mu);

      // +forward staple
      tmp1[rb[cb]] = u_nu_mu * adj(shift(u[mu],FORWARD,nu));
      tmp2[rb[cb]] = tmp1 * adj(u[nu]);

      result[rb[cb]] += param.beta * tmp2;

      // +backward staple
      tmp1[rb[cb]] = adj(shift(u_nu_mu,BACKWARD,nu)) * adj(shift(u[mu],BACKWARD,nu));
      tmp2[rb[cb]] = tmp1 * shift(u[nu],BACKWARD,nu);

      result[rb[cb]] += param.beta * tmp2;
    }

    // NOTE: a heatbath code should be responsible for resetting links on
    // a boundary. The staple is not really the correct place.

    END_CODE();
     
	    
	/*    
      QDPIO::cerr << __func__ << ": staple not possible\n";
      QDP_abort(1);*/
    }


    //! Compute dS/dU
    void RMCCoupledGaugeAct::deriv(multi1d<LatticeColorMatrix>& ds_u,
			 const Handle< GaugeState<P,Q> >& state) const
    {
      START_CODE();

      LatticeColorMatrix tmp_1;
      LatticeColorMatrix tmp_2;

      const multi1d<LatticeColorMatrix>& u = state->getLinks();
      multi1d<LatticeColorMatrix> deriv(Nd); 
      
     
      ds_u.resize(Nd);
      for(int mu=0; mu < Nd; mu++) 
      {
	deriv[mu] = zero ;

	
	for(int nu = 0; nu < Nd; nu++) 
	{ 
	  if (mu == nu) continue;

	  LatticeColorMatrix tmp_1 = shift(u[nu], FORWARD, mu);
	  LatticeColorMatrix tmp_2 = shift(u[mu], FORWARD, nu);

	  LatticeColorMatrix up_plq   = u[mu]*tmp_1*adj(tmp_2)*adj(u[nu]);
	  LatticeColorMatrix down_plq = u[mu]*shift(adj(tmp_1)*adj(u[mu])*u[nu],BACKWARD,nu);

	  deriv[mu] += up_plq + down_plq ;	                   
	                   
	
	}// nu
	// Fold in the normalization from the action
	ds_u[mu]  = (-param.beta/Real(2*Nc)        ) * deriv[mu];
		
      }// mu

      // Zero the force on any fixed boundaries
      getGaugeBC().zero(ds_u);

      END_CODE();
    }

    //! Compute the plaquette for the shadow field
    void RMCCoupledGaugeAct::siteAction_shadow(multi2d<LatticeColorMatrix>& plq) const
    {
      START_CODE();

      // Initialize
      plq.resize(Nd,Nd);
      plq = zero;

      // Handle< const GaugeState<P,Q> > u_bc(createState(u));
      // Apply boundaries
      //const multi1d<LatticeColorMatrix>& u = state->getLinks();

      // Compute the average plaquettes
      for(int mu=1; mu < Nd; ++mu)
      {
	for(int nu=0; nu < mu; ++nu)
	{
	  plq[mu][nu] += (param.shadow_gauge_field[mu]*shift(param.shadow_gauge_field[nu],FORWARD,mu)*adj(shift(param.shadow_gauge_field[mu],FORWARD,nu))*adj(param.shadow_gauge_field[nu])) ; 

	  // Keep a copy
	  plq[nu][mu] = plq[mu][nu];
	}
      }

      END_CODE();
    }


     //!  Update coupling coefficient Gamma using RMC
    void RMCCoupledGaugeAct::updateGamma(const Handle< GaugeState<P,Q> >& state) const
    {
	START_CODE();	
	//save the seed if need be
//	QDP::Seed bkup;
//	QDP::RNG::savern(bkup);
	if(((param.curr_num-1) % param.freeze_steps == 0 && param.freeze_check) || !param.freeze_check) { //if freezing is turned on and param.freeze_steps have been taken from the last update or the freezing is turned off, perform the update
	
		//construct current shadow field file name string
		std::string copy = param.shadow_gauge_id;
		copy.append("_cfg_");
		copy.append(std::to_string(param.curr_num));
		copy.append(".lime");
	
		//read shadow field into shadow_gauge_field
		XMLReader file_xml, record_xml;
		QDPFileReader to(file_xml,copy,param.file_format);
		read(to,record_xml,param.shadow_gauge_field);
		close(to);

		multi2d<LatticeColorMatrix> plq_driven;
		multi2d<LatticeColorMatrix> plq_shadow;
       		this->siteAction(plq_driven, state);
		this->siteAction_shadow(plq_shadow);

		const multi1d<LatticeColorMatrix>& u = state->getLinks();
	
		Double M = Double(0.0); //since gamma is defined to be positive, maximum of (0-gamma*A_1*A_2) must be 0.

		Double three = Nc; 

		for (int mu=1; mu<Nd; ++mu){
			for (int nu=0; nu<mu; ++nu) {
				LatticeReal plaq_driven = three-real(trace(plq_driven[mu][nu])); //calculate action density for the driven field
				LatticeReal plaq_shadow = three-real(trace(plq_shadow[mu][nu])); //calculate action density for the shadow field
			
				LatticeReal rnd_num; //fill lattice with random numbers
				random(rnd_num);
			
				//calculate RMC switching probability
				LatticeReal rnd_prob = exp((-param.gamma)/Double(Nc*Nc)*(M+plaq_shadow*plaq_driven));

				//set coupling to gamma if condition satisfied 
				param.RMCGamma[mu][nu] = where(rnd_num <= rnd_prob,param.gamma,Real(zero));
				param.RMCGamma[nu][mu] = param.RMCGamma[mu][nu]; 
			}
		}	
	}
//	QDP::RNG::setrn(bkup);
	param.curr_num++; //reflect the current update number

	END_CODE();
	
     }

     //! Restore/initialize coupling coefficient Gamma to original values
     void RMCCoupledGaugeAct::restoreGamma() const
     {
	param.RMCGamma = zero;
     }

     //! Compute the RMC Action
     Double RMCCoupledGaugeAct::RMC_S(const Handle< GaugeState<P,Q> >& state) const
     {
	START_CODE();

 	multi2d<LatticeColorMatrix> plq_driven;
	multi2d<LatticeColorMatrix> plq_shadow;
        this->siteAction(plq_driven, state);
	this->siteAction_shadow(plq_shadow);

	Double s_pg = zero;

	Double M = Double(0.0);

	Double three = Nc; 
	
 	for(int mu=1; mu < Nd; ++mu) {
		for(int nu=0; nu < mu; ++nu) {
			
			LatticeReal plaq_driven = three-real(trace(plq_driven[mu][nu]));
			LatticeReal plaq_shadow = three-real(trace(plq_shadow[mu][nu]));
			
			LatticeReal exp_term = exp((-param.gamma)/Double(Nc*Nc)*(M+plaq_shadow*plaq_driven));
		
			LatticeReal log_term = where(param.RMCGamma[mu][nu] != param.gamma, log(Real(1.0)-exp_term), Real(zero));

			LatticeReal old_s = (param.alpha/Real(Nc))*plaq_driven;

			LatticeReal gamma_s = (param.RMCGamma[mu][nu]/Double(Nc*Nc))*plaq_driven*plaq_shadow;

			s_pg += sum(old_s + gamma_s - log_term);
			
		}

	}

	END_CODE();

	return s_pg;

      }

  //! Compute dS_{RMC}/dU
  void RMCCoupledGaugeAct::RMC_deriv(multi1d<LatticeColorMatrix> & ds_u, const Handle< GaugeState<P,Q> >& state) const 
      {
	START_CODE();

    	ds_u.resize(Nd);

    	ds_u = zero;

        const multi1d<LatticeColorMatrix>& u = state->getLinks();
      
	Double M = Double(0.0);

	Double three = Nc; 

	for (int mu=0; mu<Nd; mu++) {
		for (int nu=0; nu<Nd; nu++) {

			if (mu == nu) continue;
			
			LatticeColorMatrix tmp_1 = shift(u[nu], FORWARD, mu);
	  		LatticeColorMatrix tmp_2 = shift(u[mu], FORWARD, nu);

			LatticeColorMatrix up_plaq_driven =  u[mu]*tmp_1*adj(tmp_2)*adj(u[nu]);
			LatticeColorMatrix down_plaq_driven = u[mu]*shift(adj(tmp_1)*adj(u[mu])*u[nu],BACKWARD,nu);

			LatticeColorMatrix tmp_3 = shift(param.shadow_gauge_field[nu], FORWARD, mu);
	  		LatticeColorMatrix tmp_4 = shift(param.shadow_gauge_field[mu], FORWARD, nu);

			LatticeColorMatrix up_plaq_shadow = param.shadow_gauge_field[mu]*tmp_3*adj(tmp_4)*adj(param.shadow_gauge_field[nu]);
			LatticeColorMatrix down_plaq_shadow = param.shadow_gauge_field[mu]*shift(adj(tmp_3)*adj(param.shadow_gauge_field[mu])*param.shadow_gauge_field[nu],BACKWARD,nu);

	  		LatticeColorMatrix ds_orig = -param.alpha/(Real(2*Nc))*(up_plaq_driven + down_plaq_driven);

			LatticeReal up_plaq_driven_tr = three-real(trace(up_plaq_driven));
			LatticeReal down_plaq_driven_tr = three-real(trace(down_plaq_driven));
			
			LatticeReal up_plaq_shadow_tr = three-real(trace(up_plaq_shadow));
			LatticeReal down_plaq_shadow_tr = three-real(trace(down_plaq_shadow));

			LatticeReal exp_term_up = exp((-param.gamma)/Double(Nc*Nc)*(M+up_plaq_shadow_tr*up_plaq_driven_tr));
			LatticeReal exp_term_down = exp((-param.gamma)/Double(Nc*Nc)*(M+down_plaq_shadow_tr*down_plaq_driven_tr));

			LatticeColorMatrix ds_log_up = where(param.RMCGamma[mu][nu] != param.gamma, -param.gamma*(exp_term_up/(exp_term_up-Real(1.0))) * up_plaq_shadow_tr * -up_plaq_driven/Real(2.0*Nc*Nc), ColorMatrix(zero));
			LatticeColorMatrix ds_log_down = where(shift(param.RMCGamma[mu][nu],BACKWARD,nu) != param.gamma, -param.gamma*(exp_term_down/(exp_term_down-Real(1.0))) * down_plaq_shadow_tr * -down_plaq_driven/Real(2.0*Nc*Nc), ColorMatrix(zero));

			LatticeColorMatrix ds_gamma_up = param.RMCGamma[mu][nu] * up_plaq_shadow_tr  * -up_plaq_driven/Real(2.0*Nc*Nc);

			LatticeColorMatrix ds_gamma_down = shift(param.RMCGamma[mu][nu],BACKWARD,nu) * down_plaq_shadow_tr  * -down_plaq_driven/Real(2.0*Nc*Nc);

			ds_u[mu] += ds_orig + ds_gamma_up + ds_gamma_down - ds_log_up - ds_log_down;

		}

	} 

	// Zero the force on any fixed boundaries
   	getGaugeBC().zero(ds_u);	
	
	END_CODE();

      }

  }//RMCCoupledGaugeActEnv

} // Chroma
