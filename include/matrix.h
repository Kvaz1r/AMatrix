#pragma once

#include <algorithm>
#include <array>
#include <iostream>

// A matrix Library to be simple and fast
// By Pooyan



namespace AMatrix
{



//product between two matrices
template<class T1, class T2>
class MatrixMatrixProdExpression
{
public:
    typedef typename T1::TDataType TDataType;

    MatrixMatrixProdExpression(const T1& a,
                               const T2& b):
        ma(a), mb(b)
    {}

    inline typename T1::TDataType operator()(int i, int j) const
    {
        TDataType tmp = 0.0;
        for(int k=0; k<ma.size2(); ++k)
            tmp += ma(i,k)*mb(k,j);
        return tmp;
    }

    inline int size1() const
    {
        return ma.size1();
    }
    inline int size2() const
    {
        return mb.size2();
    }

private:
    const T1& ma;
    const T2& mb;
};


//product scalar*matrix (or expression)
template<typename TScalarType, class T2>
class ScalarTimesMatrixMatrixProdExpression
{
public:
    typedef typename T2::TDataType TDataType;

    ScalarTimesMatrixMatrixProdExpression(const TScalarType s,
                                          const T2& b):
        ms(s), mb(b)
    {}

    inline typename T2::TDataType operator()(int i, int j) const
    {
        return ms*mb(i,j);
    }

    inline int size1() const
    {
        return mb.size1();
    }
    inline int size2() const
    {
        return mb.size2();
    }

private:
    const TScalarType ms;
    const T2& mb;
};

//utility to set to zero the matrix
template<typename TData, int TSize1, int TSize2>
class ZeroMatrixExpression
{
public:
    typedef TData TDataType;

    ZeroMatrixExpression()
    {}

    inline TDataType operator()(int i, int j) const
    {
        return 0.0;
    }

    inline constexpr int size1() const
    {
        return TSize1;
    }
    inline constexpr int size2() const
    {
        return TSize2;
    }

};


template <typename DataType, int NumberOfRows, int NumberOfColumns>
class Matrix
{

    DataType _data[NumberOfRows * NumberOfColumns];

public:
    typedef DataType TDataType;

    Matrix() {}
    Matrix(DataType const& InitialValue)
    {
        for (int i = 0; i < size(); i++)
            _data[i] = InitialValue;
    }

    Matrix(Matrix const& Other)
    {
        for (int i = 0; i < size(); i++)
            _data[i] = Other._data[i];
    }

    Matrix(Matrix&& Other) = default;

    template <typename TOtherMatrixType>
    Matrix(TOtherMatrixType const& Other)
    {
        for (int i = 0; i < size1(); i++)
            for (int j = 0; j < size2(); j++)
                at(i, j) = Other(i, j);
    }

    Matrix& operator=(Matrix const& Other)
    {
        for (int i = 0; i < size(); i++)
            _data[i] = Other._data[i];
        return *this;
    }

    template<class T>
    Matrix& operator=(T const& Other)
    {
        for(int i=0; i<size1(); ++i)
            for(int j=0; j<size2(); ++j)
                (*this)(i,j) = Other(i,j);
        return *this;
    }

    DataType& operator()(int i, int j)
    {
        return at(i, j);
    }

    DataType const& operator()(int i, int j) const
    {
        return at(i, j);
    }

    DataType& at(int i, int j)
    {
        return _data[i * NumberOfColumns + j];
    }

    DataType const& at(int i, int j) const
    {
        return _data[i * NumberOfColumns + j];
    }

    constexpr int size1() const
    {
        return NumberOfRows;
    }

    constexpr int size2() const
    {
        return NumberOfColumns;
    }

    constexpr int size() const
    {
        return NumberOfRows * NumberOfColumns;
    }

    friend bool operator==(Matrix const& First, Matrix const& Second)
    {
        for (int i = 0; i < First.size(); i++)
            if (First._data[i] != Second._data[i])
                return false;
        return true;
    }

    friend Matrix operator+(Matrix const& First, Matrix const& Second)
    {
        Matrix result;
        const DataType* __restrict first_data = First._data;
        const DataType* __restrict second_data = Second._data;
        DataType* __restrict result_data = result._data;
        for (int i = 0; i < First.size(); ++i)
            *result_data++ = *first_data++ + *second_data++;

        return result;
    }

    friend Matrix operator-(Matrix const& First, Matrix const& Second)
    {
        Matrix result;
        for (int i = 0; i < First.size(); ++i)
            result._data[i] = First._data[i] - Second._data[i];

        return result;
    }
    
    template< class T >
    Matrix& operator+=(const T& Other)
    {
        for(int i=0; i<size1(); ++i)
            for(int j=0; j<size2(); ++j)
                at(i,j) += Other(i,j);
            
        return *this;
    }

//     template <int SecondNumberOfColumns>
//     friend inline Matrix<DataType, NumberOfRows, SecondNumberOfColumns>
//     operator*(Matrix const& First,
//         Matrix<DataType, NumberOfColumns, SecondNumberOfColumns> const&
//             Second) {
//         Matrix<DataType, NumberOfRows, SecondNumberOfColumns> result;
//         for (int i = 0; i < NumberOfRows; i++)
//             for (int j = 0; j < SecondNumberOfColumns; j++) {
//                 DataType temp = DataType();
//                 for (int k = 0; k < NumberOfColumns; k++)
//                     temp += First(i, k) * Second(k, j);
//
//                 result(i, j) = temp;
//             }
//
//         return result;
//     }

    /*    template <class T1, class T2>
        inline MatrixMatrixProdExpression<T1,T2> operator*(const T1& First, const T2& Second) {
            return MatrixMatrixProdExpression<T1,T2>(First,Second);
        }*/

    Matrix& noalias()
    {
        return *this;
    }

    DataType* data()
    {
        return _data;
    }

    DataType const* data() const
    {
        return _data;
    }

private:
    template <int TSize>
    inline static void ElementwiseMult(const DataType* __restrict A,
                                       const DataType* __restrict B, DataType* C)
    {
        for (int i = 0; i < TSize; ++i)
        {
            *(C++) += *(A++) * *(B++);
        }
    }
};

template <class T1, class T2>
inline MatrixMatrixProdExpression<T1,T2> prod(const T1& First, const T2& Second)
{
    return MatrixMatrixProdExpression<T1,T2>(First,Second);
}

template <class T1, class T2>
inline MatrixMatrixProdExpression<T1,T2> operator*(const T1& First, const T2& Second)
{
    return MatrixMatrixProdExpression<T1,T2>(First,Second);
}

template <class T2>
inline ScalarTimesMatrixMatrixProdExpression<double,T2> operator*(const double Scalar, const T2& Second)
{
    return ScalarTimesMatrixMatrixProdExpression<double,T2>(Scalar,Second);
}

template <typename DataType, int NumberOfRows, int NumberOfColumns>
bool operator!=(Matrix<DataType, NumberOfRows, NumberOfColumns> const& First,
                Matrix<DataType, NumberOfRows, NumberOfColumns> const& Second)
{
    return !(First == Second);
}

template <class T1, class T2>
void operator+=(const T2& First, const T2& Second)
{
    for(int i=0; i<First.size1(); ++i)
        for(int j=0; j<First.size2(); ++j)
            First(i,j) += Second(i,j);
}

/// output stream function
template <typename DataType, int NumberOfRows, int NumberOfColumns>
inline std::ostream& operator<<(std::ostream& rOStream,
                                Matrix<DataType, NumberOfRows, NumberOfColumns> const& TheMatrix)
{
    rOStream << "{";
    for (int i = 0; i < NumberOfRows; i++)
    {
        for (int j = 0; j < NumberOfColumns; j++)
            rOStream << TheMatrix(i, j) << ",";
        rOStream << std::endl;
    }
    rOStream << "}";

    return rOStream;
}
}  // mamespace AMatrix
