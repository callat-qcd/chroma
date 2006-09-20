#include "prec_logdet_ee_monomial_w.h"

#include "update/molecdyn/monomial/monomial_factory.h"

#include "actions/ferm/fermacts/fermact_factory_w.h"
#include "actions/ferm/fermacts/fermacts_aggregate_w.h"

//#include "actions/ferm/fermacts/prec_clover_fermact_w.h"

using namespace std;

namespace Chroma 
{ 

  namespace PrecLogDetEvenEvenMonomial4DEnv 
  {
    namespace
    {
      //! Callback function for the factory
      Monomial< multi1d<LatticeColorMatrix>,
		multi1d<LatticeColorMatrix> >* createMonomial(XMLReader& xml, const string& path) 
      {
	QDPIO::cout << "Create Monomial: " << name << endl;

	return new PrecLogDetEvenEvenMonomial4D(PrecLogDetEvenEvenMonomialParams(xml, path));
      }

      //! Local registration flag
      bool registered = false;
    }

    const std::string name = std::string("N_FLAVOR_LOGDET_EVEN_EVEN_FERM_MONOMIAL");
 
    //! Register all the factories
    bool registerAll() 
    {
      bool success = true; 
      if (! registered)
      {
	success &= TheMonomialFactory::Instance().registerObject(name, createMonomial);
	registered = true;
      }
      return success;
    }
  }
 

 
  PrecLogDetEvenEvenMonomialParams::PrecLogDetEvenEvenMonomialParams(XMLReader& in, const std::string& path) 
  {
    XMLReader paramtop(in, path);
    
    fermact = readXMLGroup(paramtop, "FermionAction", "FermAct");
    read(paramtop,"num_flavors", num_flavors);

    QDPIO::cout << "PrecLogDetEvenEvenMonomialParams: read \n" << fermact.id << endl;
  }

  void read(XMLReader& r, const std::string& path,  PrecLogDetEvenEvenMonomialParams& p) 
  {
    PrecLogDetEvenEvenMonomialParams tmp(r, path);
    p = tmp;
  }

  void write(XMLWriter& xml, const std::string& path, const PrecLogDetEvenEvenMonomialParams& p)
  {
    // Not implemented
  }


  PrecLogDetEvenEvenMonomial4D::PrecLogDetEvenEvenMonomial4D(const PrecLogDetEvenEvenMonomialParams& p) : 
    num_flavors(p.num_flavors) 
  {
    START_CODE();

    // Grok the fermact out of the XML
    std::istringstream is(p.fermact.xml);
    XMLReader fermact_reader(is);
    QDPIO::cout << "EvanOddPrecConstDetTwoFlavorWilsonTypeFermMonomial: construct " << p.fermact.id << endl;

    WilsonTypeFermAct<T,P,Q>* tmp_act = 
      TheWilsonTypeFermActFactory::Instance().createObject(p.fermact.id, fermact_reader, p.fermact.path);

    EvenOddPrecLogDetWilsonTypeFermAct<T,P,Q>* downcast = 
      dynamic_cast<EvenOddPrecLogDetWilsonTypeFermAct<T,P,Q>*>(tmp_act);

    // Check success of the downcast 
    if( downcast == 0x0 ) {
      QDPIO::cerr << "Unable to downcast FermAct to EvenOddPrecWilsonTypeFermAct in EvenOddPrecConstDetTwoFlavorWilsonTypeFermMonomial()" << endl;
      QDP_abort(1);
    }
    
    fermact = downcast; 
    
    END_CODE();
  }


} // End namespace
