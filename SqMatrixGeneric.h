/*
Copyright 2011, Jernej Kovacic

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

   http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
*/


/**
@file SqMatrixGeneric.h

Declaration of the class SqMatrixGeneric, representing square matrices.
The class is derived from MatrixGeneric and has additional functions defined
for square matrices only

@author Jernej Kovacic
*/

#ifndef _SQMATRIXGENERIC_H_
#define	_SQMATRIXGENERIC_H_

#include "MatrixGeneric.h"

template<class T>
class SqMatrixGeneric : public MatrixGeneric<T>
{
public:
    // Constructor
    SqMatrixGeneric(unsigned int dim = 1) throw(MatrixException);
    // Copy constructor
    SqMatrixGeneric(const MatrixGeneric<T>& orig) throw(MatrixException);
    // operator= (it must be reimplemnted as it is not inherited from the base class)
    SqMatrixGeneric<T>& operator= (const MatrixGeneric<T>& m) throw (MatrixException);
    
    // Several methods to create diagonal matrices
    MatrixGeneric<T>& setDiag(const T& scalar) throw(MatrixException);
    MatrixGeneric<T>& setUnit() throw(MatrixException);

    // Determinant of the matrix
    T determinant() const throw(MatrixException);

    // Inverse matrix
    SqMatrixGeneric<T> inverse() const throw(MatrixException);

    // The following functions inherited from the parent class
    // are not applicable for square matrices and throw an exception
    MatrixGeneric<T>& removeRow(unsigned int rowNr) throw (MatrixException);
    MatrixGeneric<T>& removeColumn(unsigned int colNr) throw (MatrixException);
    MatrixGeneric<T>& insertRow(unsigned int rowNr) throw (MatrixException);
    MatrixGeneric<T>& insertColumn(unsigned int colNr) throw (MatrixException);
};

// Matrices with elements of types float and double make most sense
// so these two types are predefined
typedef SqMatrixGeneric<float> FSqMatrix;
typedef SqMatrixGeneric<double> SqMatrix;

// DEFINITION

// This is a templated class, so its definition must follow its declaration.
// When building, THIS file must be compiled.
// Alternatively the definition can be included into this file.
#include "SqMatrixGeneric.cpp"

#endif	/* _SQMATRIXGENERIC_H_ */
