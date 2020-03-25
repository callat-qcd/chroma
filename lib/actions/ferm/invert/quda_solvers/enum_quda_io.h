#ifndef enum_quda_io_h
#define enum_quda_io_h

#include "chromabase.h"
#include <string>
#include "singleton.h"
#include "io/enum_io/enum_type_map.h"

namespace Chroma
{
    /*!
     * Types and structures
     *
     * \ingroup io
     *
     * @{
     */

    //! Quda Solver type
    enum QudaSolverType { 
        CG,
        BICGSTAB,
        GCR,
        MR
    };

    namespace QudaSolverTypeEnv { 
        extern const std::string typeIDString;
        extern bool registered; 
        bool registerAll(void);   // Forward declaration
    }

    // A singleton to hold the typemap
    typedef SingletonHolder<EnumTypeMap<QudaSolverType> > theQudaSolverTypeMap;

    // Reader and writer

    //! Read an QudaSolverType enum
    void read(XMLReader& r, const std::string& path, QudaSolverType& t);

    //! Write an QudaSolverType enum
    void write(XMLWriter& w, const std::string& path, const QudaSolverType& t);


    // QUDA Lanczos type and ENV
    /** S=smallest L=largest                                                                                                                                                        
        R=real M=modulus I=imaniary **/
    enum ChromaQudaEigSpectrumType {
        SR,
        LR,
        SM,
        LM,
        SI,
        LI
    };

    namespace ChromaQudaEigSpectrumTypeEnv { 
        extern const std::string typeIDString;
        extern bool registered; 
        bool registerAll(void);   // Forward declaration
    }

    // A singleton to hold the typemap
    typedef SingletonHolder<EnumTypeMap<ChromaQudaEigSpectrumType> > theChromaQudaEigSpectrumTypeMap;

    //! Read an QudaEigSpectrumType enum
    void read(XMLReader& r, const std::string& path, ChromaQudaEigSpectrumType& t);

    //! Write an QudaEigSpectrumType enum
    void write(XMLWriter& w, const std::string& path, const ChromaQudaEigSpectrumType& t);

  

    //! Quda Precision type
    enum QudaPrecisionType { 
        DEFAULT,
        HALF,
        SINGLE,
        DOUBLE
    };
  
    namespace QudaPrecisionTypeEnv { 
        extern const std::string typeIDString;
        extern bool registered; 
        bool registerAll(void);   // Forward declaration
    }

    // A singleton to hold the typemap
    typedef SingletonHolder<EnumTypeMap<QudaPrecisionType> > theQudaPrecisionTypeMap;

    // Reader and writer

    //! Read an QudaSolverType enum
    void read(XMLReader& r, const std::string& path, QudaPrecisionType& t);

    //! Write an QudaSolverType enum
    void write(XMLWriter& w, const std::string& path, const QudaPrecisionType& t);


    //! Quda Gauge Reconstruct type
    enum QudaReconsType { 
        RECONS_NONE,
        RECONS_8,
        RECONS_12
    };
  
    namespace QudaReconsTypeEnv { 
        extern const std::string typeIDString;
        extern bool registered; 
        bool registerAll(void);   // Forward declaration
    }

    // A singleton to hold the typemap
    typedef SingletonHolder<EnumTypeMap<QudaReconsType> > theQudaReconsTypeMap;

    // Reader and writer

    //! Read an QudaReconsType enum
    void read(XMLReader& r, const std::string& path, QudaReconsType& t);

    //! Write an QudaReconsType enum
    void write(XMLWriter& w, const std::string& path, const QudaReconsType& t);


    enum QudaSchwarzMethod { 
        ADDITIVE_SCHWARZ,
        MULTIPLICATIVE_SCHWARZ
    };
  
    namespace QudaSchwarzMethodEnv { 
        extern const std::string typeIDString;
        extern bool registered; 
        bool registerAll(void);   // Forward declaration
    }

    // A singleton to hold the typemap
    typedef SingletonHolder<EnumTypeMap<QudaSchwarzMethod> > theQudaSchwarzMethodMap;

    // Reader and writer

    //! Read an QudaSchwarzMethod enum
    void read(XMLReader& r, const std::string& path, QudaSchwarzMethod& t);

    //! Write an QudaSchwarzMethod enum
    void write(XMLWriter& w, const std::string& path, const QudaSchwarzMethod& t);



    /*! @} */   // end of group io






}




#endif
