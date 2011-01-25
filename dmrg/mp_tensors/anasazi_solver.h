#ifndef ANASAZI_MULTI_VEC_TRAITS_HPP
#define ANASAZI_MULTI_VEC_TRAITS_HPP

#include "AnasaziConfigDefs.hpp"
#include "AnasaziTypes.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_SerialDenseMatrix.hpp"

#include <boost/shared_ptr.hpp>
#include <boost/noncopyable.hpp>

namespace Anasazi {
    
    // They want a shallow copy...what a bunch of morons
    // To avoid confusion, I derive from noncopyable
    // To implement it, I use a shared_ptr
    
    template<class Vector>
    class IETLMultMv : public boost::noncopyable
    {
    public:
        typedef std::vector<Vector> data_type;
        
        IETLMultMv(std::vector<Vector> const & v)
        : data( new std::vector<Vector>(v) )
        { }
        
        IETLMultMv() { }
        
        boost::shared_ptr<std::vector<Vector> > data;
    }
    
    template< class ScalarType, class MV >
    struct UndefinedMultiVecTraits
    {
        static inline ScalarType notDefined() { return MV::this_type_is_missing_a_specialization(); };
    };
    
    template<class ScalarType, class Vector>
    class MultiVecTraits<double, IETLMultMv<Vector> >
    {
        typedef IETLMultMv<Vector> MV;
        
    public:
        static Teuchos::RCP<MV> Clone( const MV& mv, const int numvecs )
        {
            return Teuchos::RCP<MV>(new MV( typename MV::data_type(numvecs) ));
        }
        
        static Teuchos::RCP<MV> CloneCopy( const MV& mv )
        {
            return Teuchos::RCP<MV>(new MV( mv.data ));
        }
        
        static Teuchos::RCP<MV> CloneCopy( const MV& mv, const std::vector<int>& index )
        {
            MV * ret = new MV;
            for (std::vector<int>::const_iterator it = index.begin();
                 it != index.end(); ++it)
                ret->data->push_back( (*mv.data)[*it] );
            return Teuchos::RCP<MV>(ret);
        }
        
        static Teuchos::RCP<MV> CloneViewNonConst( MV& mv, const std::vector<int>& index )
        {
            MV *
        
        static Teuchos::RCP<const MV> CloneView( const MV& mv, const std::vector<int>& index )
        { UndefinedMultiVecTraits<ScalarType, MV>::notDefined(); return Teuchos::null; }
        
        
        
        static int GetVecLength( const MV& mv )
        { UndefinedMultiVecTraits<ScalarType, MV>::notDefined(); return 0; }     
        
        static int GetNumberVecs( const MV& mv )
        { UndefinedMultiVecTraits<ScalarType, MV>::notDefined(); return 0; }     
        
        
        
        
        static void MvTimesMatAddMv( const ScalarType alpha, const MV& A, 
                                    const Teuchos::SerialDenseMatrix<int,ScalarType>& B, 
                                    const ScalarType beta, MV& mv )
        { UndefinedMultiVecTraits<ScalarType, MV>::notDefined(); }     
        
        static void MvAddMv( const ScalarType alpha, const MV& A, const ScalarType beta, const MV& B, MV& mv )
        { UndefinedMultiVecTraits<ScalarType, MV>::notDefined(); }     
        
        static void MvScale ( MV& mv, const ScalarType alpha )
        { UndefinedMultiVecTraits<ScalarType, MV>::notDefined(); }
        
        static void MvScale ( MV& mv, const std::vector<ScalarType>& alpha )
        { UndefinedMultiVecTraits<ScalarType, MV>::notDefined(); }
        
        static void MvTransMv( const ScalarType alpha, const MV& A, const MV& mv, Teuchos::SerialDenseMatrix<int,ScalarType>& B)
        { UndefinedMultiVecTraits<ScalarType, MV>::notDefined(); }     
        
        static void MvDot ( const MV& mv, const MV& A, std::vector<ScalarType> &b) 
        { UndefinedMultiVecTraits<ScalarType, MV>::notDefined(); }     
        
        
        
        
        static void MvNorm( const MV& mv, std::vector<typename Teuchos::ScalarTraits<ScalarType>::magnitudeType> &normvec )
        { UndefinedMultiVecTraits<ScalarType, MV>::notDefined(); }     
        
        
        
        
        static void SetBlock( const MV& A, const std::vector<int>& index, MV& mv )
        { UndefinedMultiVecTraits<ScalarType, MV>::notDefined(); }     
        
        static void MvRandom( MV& mv )
        { UndefinedMultiVecTraits<ScalarType, MV>::notDefined(); }     
        
        static void MvInit( MV& mv, const ScalarType alpha = Teuchos::ScalarTraits<ScalarType>::zero() )
        { UndefinedMultiVecTraits<ScalarType, MV>::notDefined(); }     
        
        
        
        
        static void MvPrint( const MV& mv, std::ostream& os )
        { UndefinedMultiVecTraits<ScalarType, MV>::notDefined(); }     
        
    };
    
} // namespace Anasazi

#endif // ANASAZI_MULTI_VEC_TRAITS_HPP