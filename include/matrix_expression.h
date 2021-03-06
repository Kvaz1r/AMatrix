#pragma once

#include <iostream>

namespace AMatrix {
constexpr std::size_t dynamic = 0;
constexpr std::size_t row_major_access = 1;
constexpr std::size_t column_major_access = 2;
constexpr std::size_t unordered_access = 3;

template <std::size_t TCategory1, std::size_t TCategory2>
class AccessTrait {
   public:
    static constexpr std::size_t Category = unordered_access;
};

template <>
class AccessTrait<row_major_access, row_major_access> {
   public:
    static constexpr std::size_t category = row_major_access;
};

template <typename TExpressionType, std::size_t TCategory = unordered_access>
class MatrixExpression {
   public:
    static constexpr std::size_t category = TCategory;

    // using value_type = TExpressionType::value_type;
    MatrixExpression() {}

    MatrixExpression(MatrixExpression const& Other) = default;

    MatrixExpression(MatrixExpression&& Other) = default;

    MatrixExpression& operator=(MatrixExpression const& Other) = default;

    MatrixExpression& operator=(MatrixExpression&& Other) = default;

    TExpressionType& expression() {
        return *static_cast<TExpressionType*>(this);
    }
    TExpressionType const& expression() const {
        return *static_cast<const TExpressionType* const>(this);
    }

    TExpressionType& noalias() { return expression(); }
};

template <typename TExpressionType>
class TransposeMatrix
    : public MatrixExpression<TransposeMatrix<TExpressionType>> {
    TExpressionType const& mOriginal;

   public:
    using data_type = typename TExpressionType::data_type;
    TransposeMatrix() = delete;

    TransposeMatrix(TExpressionType const& Original) : mOriginal(Original) {}

    inline data_type const& operator()(std::size_t i, std::size_t j) const {
        return mOriginal(j, i);
    }

    inline std::size_t size1() const { return mOriginal.size2(); }
    inline std::size_t size2() const { return mOriginal.size1(); }
};

template <typename TExpressionType>
class SubMatrix : public MatrixExpression<SubMatrix<TExpressionType>> {
    TExpressionType& mOriginal;
    std::size_t _origin_index1;
    std::size_t _origin_index2;
    std::size_t _size1;
    std::size_t _size2;

   public:
    using data_type = typename TExpressionType::data_type;
    SubMatrix() = delete;

    SubMatrix(TExpressionType& Original, std::size_t OriginIndex1,
        std::size_t OriginIndex2, std::size_t TheSize1, std::size_t TheSize2)
        : mOriginal(Original),
          _origin_index1(OriginIndex1),
          _origin_index2(OriginIndex2),
          _size1(TheSize1),
          _size2(TheSize2) {}

    inline data_type const& operator()(std::size_t i, std::size_t j) const {
        return mOriginal(i + _origin_index1, j + _origin_index2);
    }

    inline data_type& operator()(std::size_t i, std::size_t j) {
        return mOriginal(i + _origin_index1, j + _origin_index2);
    }

    inline std::size_t size1() const { return _size1; }
    inline std::size_t size2() const { return _size2; }
};

template <typename TDataType>
class ZeroMatrix
    : public MatrixExpression<ZeroMatrix<TDataType>, row_major_access> {
    std::size_t _size1;
    std::size_t _size2;

   public:
    using data_type = TDataType;

    ZeroMatrix() = delete;

    ZeroMatrix(std::size_t Size) : _size1(Size), _size2(1) {}

    ZeroMatrix(std::size_t Size1, std::size_t Size2)
        : _size1(Size1), _size2(Size2) {}

    inline TDataType operator()(std::size_t i, std::size_t j) const {
        return TDataType();
    }

    inline TDataType operator[](std::size_t i) const { return TDataType(); }

    inline std::size_t size1() const { return _size1; }
    inline std::size_t size2() const { return _size2; }
};

template <typename TExpression1Type, typename TExpression2Type>
class MatrixSumExpression
    : public MatrixExpression<
          MatrixSumExpression<TExpression1Type, TExpression2Type>,
          AccessTrait<TExpression1Type::category,
              TExpression2Type::category>::category> {
    TExpression1Type const& _first;
    TExpression2Type const& _second;

   public:
    MatrixSumExpression(
        TExpression1Type const& First, TExpression2Type const& Second)
        : _first(First), _second(Second) {}
    using data_type = typename TExpression1Type::data_type;

    std::size_t size1() const { return _first.size1(); }

    std::size_t size2() const { return _first.size2(); }

    std::size_t size() const { return _first.size(); }

    inline data_type operator()(std::size_t i, std::size_t j) const {
        return _first(i, j) + _second(i, j);
    }

    inline data_type operator[](std::size_t i) const {
        return _first[i] + _second[i];
    }
};

template <typename TExpression1Type, typename TExpression2Type,
    std::size_t TCategory1, std::size_t TCategory2>
MatrixSumExpression<TExpression1Type, TExpression2Type> operator+(
    MatrixExpression<TExpression1Type, TCategory1> const& First,
    MatrixExpression<TExpression2Type, TCategory2> const& Second) {
    return MatrixSumExpression<TExpression1Type, TExpression2Type>(
        First.expression(), Second.expression());
}

template <typename TExpression1Type, typename TExpression2Type>
class MatrixMinusExpression
    : public MatrixExpression<
          MatrixMinusExpression<TExpression1Type, TExpression2Type>,
          AccessTrait<TExpression1Type::category,
              TExpression2Type::category>::category> {
    TExpression1Type const& _first;
    TExpression2Type const& _second;

   public:
    MatrixMinusExpression(
        TExpression1Type const& First, TExpression2Type const& Second)
        : _first(First), _second(Second) {}
    using data_type = typename TExpression1Type::data_type;

    std::size_t size1() const { return _first.size1(); }

    std::size_t size2() const { return _first.size2(); }

    std::size_t size() const { return _first.size(); }

    inline data_type operator()(std::size_t i, std::size_t j) const {
        return _first(i, j) - _second(i, j);
    }

    inline data_type operator[](std::size_t i) const {
        return _first[i] - _second[i];
    }
};

template <typename TExpression1Type, typename TExpression2Type,
    std::size_t TCategory1, std::size_t TCategory2>
MatrixMinusExpression<TExpression1Type, TExpression2Type> operator-(
    MatrixExpression<TExpression1Type, TCategory1> const& First,
    MatrixExpression<TExpression2Type, TCategory2> const& Second) {
    return MatrixMinusExpression<TExpression1Type, TExpression2Type>(
        First.expression(), Second.expression());
}

template <typename TExpressionType>
class MatrixScalarProductExpression
    : public MatrixExpression<MatrixScalarProductExpression<TExpressionType>,
          AccessTrait<TExpressionType::category, row_major_access>::category> {
    typename TExpressionType::data_type const& _first;
    TExpressionType const& _second;

   public:
    using data_type = typename TExpressionType::data_type;

    MatrixScalarProductExpression(
        data_type const& First, TExpressionType const& Second)
        : _first(First), _second(Second) {}
    std::size_t size1() const { return _second.size1(); }

    std::size_t size2() const { return _second.size2(); }

    std::size_t size() const { return _second.size(); }

    inline data_type operator()(std::size_t i, std::size_t j) const {
        return _first * _second(i, j);
    }

    inline data_type operator[](std::size_t i) const {
        return _first * _second[i];
    }
};

template <typename TExpressionType, std::size_t TCategory>
MatrixScalarProductExpression<TExpressionType> operator*(
    typename TExpressionType::data_type const& First,
    MatrixExpression<TExpressionType, TCategory> const& Second) {
    return MatrixScalarProductExpression<TExpressionType>(
        First, Second.expression());
}

template <typename TExpressionType, std::size_t TCategory>
MatrixScalarProductExpression<TExpressionType> operator*(
    MatrixExpression<TExpressionType, TCategory> const& First,
    typename TExpressionType::data_type const& Second) {
    return MatrixScalarProductExpression<TExpressionType>(
        Second, First.expression());
}

template <typename TExpression1Type, typename TExpression2Type>
class MatrixProductExpression
    : public MatrixExpression<
          MatrixProductExpression<TExpression1Type, TExpression2Type>,
          unordered_access> {
    TExpression1Type const& _first;
    TExpression2Type const& _second;

   public:
    MatrixProductExpression(
        TExpression1Type const& First, TExpression2Type const& Second)
        : _first(First), _second(Second) {}
    using data_type = typename TExpression1Type::data_type;

    std::size_t size1() const { return _first.size1(); }

    std::size_t size2() const { return _second.size2(); }

    std::size_t size() const { return size1() * size2(); }

    inline data_type operator()(std::size_t i, std::size_t j) const {
        data_type result = data_type();
        for (std::size_t k = 0; k < _first.size2(); k++)
            result += _first(i, k) * _second(k, j);
        return result;
    }
};

template <typename TExpression1Type, typename TExpression2Type,
    std::size_t TCategory1, std::size_t TCategory2>
MatrixProductExpression<TExpression1Type, TExpression2Type> operator*(
    MatrixExpression<TExpression1Type, TCategory1> const& First,
    MatrixExpression<TExpression2Type, TCategory2> const& Second) {
    return MatrixProductExpression<TExpression1Type, TExpression2Type>(
        First.expression(), Second.expression());
}

}  // namespace AMatrix