/*
 * This file is part of the repicea-statistics library.
 *
 * Copyright (C) 2009-2014 Mathieu Fortin for Rouge-Epicea
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 3 of the License, or (at your option) any later version.
 *
 * This library is distributed with the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied
 * warranty of MERCHANTABILITY or FITNESS FOR A
 * PARTICULAR PURPOSE. See the GNU Lesser General Public
 * License for more details.
 *
 * Please see the license at http://www.gnu.org/copyleft/lesser.html.
 */
package repicea.stats.distributions;

import java.security.InvalidParameterException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import repicea.math.DiagonalMatrix;
import repicea.math.Matrix;
import repicea.math.SymmetricMatrix;
import repicea.stats.Distribution;
import repicea.stats.StatisticalUtility;
import repicea.stats.StatisticalUtility.TypeMatrixR;

/**
 * The CenteredGaussianDistribution class is designed for dealing with univariate or multivariate 
 * error terms in linear and non linear model. It assumes a constant variance that is represented 
 * by the variance parameter in the constructor.
 * @author Mathieu Fortin - August 2014
 */
@SuppressWarnings("serial")
public final class CenteredGaussianDistribution implements Distribution<Matrix, SymmetricMatrix> {
	
	private final GaussianDistribution underlyingDistribution;
	private final double correlationParameter;
	private final TypeMatrixR type;
	private final boolean isStructured;
	
	private final Map<List<Integer>, SymmetricMatrix> structuredVarianceCovarianceMap;
	private final Map<List<Integer>, Matrix> structuredLowerCholeskyMap;
	private final Map<Integer, SymmetricMatrix> simpleVarianceCovarianceMap;
	private final Map<Integer, Matrix> simpleLowerCholeskyMap;
	
	/**
	 * General constructor.
	 * @param variance the homogeneous variance
	 * @param correlationParameter the correlation parameter in the correlation structure
	 * @param type a TypeMatrixR enum
	 */
	public CenteredGaussianDistribution(SymmetricMatrix variance, double correlationParameter, TypeMatrixR type) {
		underlyingDistribution = new GaussianDistribution(new Matrix(variance.m_iRows,1), variance);
		this.correlationParameter = correlationParameter;
		this.type = type;
		if (this.correlationParameter != 0 && this.type != null) {
			isStructured = true;
			if (variance.m_iRows > 1) {
				throw new InvalidParameterException("The CenteredGaussianDistribution is not designed for a multivariate distribution with heterogeneous variances yet.");
			}
		} else {
			isStructured = false;
		}
		structuredVarianceCovarianceMap = new HashMap<List<Integer>, SymmetricMatrix>();
		structuredLowerCholeskyMap = new HashMap<List<Integer>, Matrix>();
		simpleVarianceCovarianceMap = new HashMap<Integer, SymmetricMatrix>();
		simpleLowerCholeskyMap = new HashMap<Integer, Matrix>();
	}

	/**
	 * Constructor for univariate distribution.
	 * @param variance the homogeneous variance
	 */
	public CenteredGaussianDistribution(SymmetricMatrix variance) {
		this(variance, 0, null);
	}

	private Matrix getLowerCholesky(List<Integer> indexList) {
		if (isStructured()) {
			if (!structuredLowerCholeskyMap.containsKey(indexList)) {
				updateMaps(indexList);
			}
			return structuredLowerCholeskyMap.get(indexList);
		} else {
			int size = indexList.size();
			if (!simpleLowerCholeskyMap.containsKey(size)) {
				updateMaps(size);
			}
			return simpleLowerCholeskyMap.get(size);
		}
	}
		
	@SuppressWarnings({ "unchecked", "rawtypes" })
	private void updateMaps(Object key) {
		if (key instanceof List) {
			List<Integer> referenceList = new ArrayList<Integer>();
			referenceList.addAll(((List) key));		// make a copy to avoid changes through reference
			Matrix distances = new Matrix(referenceList);
//			Matrix correlationMatrix = StatisticalUtility.constructRMatrix(distances, 1d, correlationParameter, type);
			SymmetricMatrix correlationMatrix = StatisticalUtility.constructRMatrix(Arrays.asList(new Double[] {1d, correlationParameter}), type, distances);
			SymmetricMatrix varianceCovariance = correlationMatrix.scalarMultiply(underlyingDistribution.getVariance().getValueAt(0, 0));
			structuredVarianceCovarianceMap.put(referenceList, varianceCovariance);
			Matrix lowerChol = varianceCovariance.getLowerCholTriangle();
			structuredLowerCholeskyMap.put(referenceList, lowerChol);
		} else {
			int size = (Integer) key;
			DiagonalMatrix varianceCovariance = Matrix.getIdentityMatrix(size).scalarMultiply(underlyingDistribution.getVariance().getValueAt(0, 0));
			simpleVarianceCovarianceMap.put(size, varianceCovariance);
			Matrix lowerChol = varianceCovariance.getLowerCholTriangle();
			simpleLowerCholeskyMap.put(size, lowerChol);
		}
	}

	private SymmetricMatrix getVariance(List<Integer> indexList) {
		if (isStructured()) {
			if (!structuredVarianceCovarianceMap.containsKey(indexList)) {
				updateMaps(indexList);
			}
			return structuredVarianceCovarianceMap.get(indexList);
		} else {
			int size = indexList.size();
			if (!simpleVarianceCovarianceMap.containsKey(size)) {
				updateMaps(size);
			}
			return simpleVarianceCovarianceMap.get(size);
		}
	}
	
	@Override
	public Matrix getMean() {
		return underlyingDistribution.getMean();
	}

	/**
	 * This method should be used in preference to the getMean() method.
	 * @param errorTermList a GaussianErrorTermList instance
	 * @return a Matrix instance
	 */
	public Matrix getMean(GaussianErrorTermList errorTermList) {
		if (errorTermList != null & !errorTermList.isEmpty()) {
			Matrix chol = getLowerCholesky(errorTermList.getDistanceIndex());
			return chol.multiply(errorTermList.getNormalizedErrors());
		} else {
			return null;
		}
	}
	
	@Override
	public SymmetricMatrix getVariance() {
		return underlyingDistribution.getVariance();
	}

	/**
	 * Provide the variance of the distribution given some error terms.<p>
	 * The class adapts the variance matrix as the number of error terms increases.
	 * @param errorTermList a GaussianErrorTermList instance
	 * @return a SymmetricMatrix instance
	 */
	public SymmetricMatrix getVariance(GaussianErrorTermList errorTermList) {
		if (errorTermList != null & !errorTermList.isEmpty()) {
			return getVariance(errorTermList.getDistanceIndex());
		} else {
			return null;
		}
	}
	
	@Override
	public Matrix getRandomRealization() {
		return underlyingDistribution.getRandomRealization();
	}

	public Matrix getRandomRealization(GaussianErrorTermList errorTermList) {
		if (errorTermList != null & !errorTermList.isEmpty()) {
			Matrix errorTerms;
			if (!errorTermList.updated) {
				List<Integer> indexList = errorTermList.getDistanceIndex();
				Matrix stdMat = getLowerCholesky(indexList);
				Matrix normalizedErrorTerms = errorTermList.getNormalizedErrors();
				errorTerms = stdMat.multiply(normalizedErrorTerms);
				errorTermList.updateErrorTerm(errorTerms);
			} else {
				errorTerms = errorTermList.getRealizedErrors();
			}
			return errorTerms;
		} else {
			return null;
		}
	}
	
	
	
	@Override
	public boolean isParametric() {return true;}

	@Override
	public boolean isMultivariate() {return getMean().m_iRows > 1;}

	@Override
	public Type getType() {return Type.GAUSSIAN;}
	
	public boolean isStructured() {return isStructured;}
	
//	public static void main(String[] args) {
//		List<Integer> list1 = new ArrayList<Integer>();
//		list1.add(1);
//		list1.add(2);
//		List<Integer> list2 = new ArrayList<Integer>();
//		list2.add(2);
//		list2.add(1);
//		System.out.println("Lists are equal : " + list1.equals(list2));
//		Collections.sort(list2);
//		System.out.println("Lists are equal : " + list1.equals(list2));
//	
//	}

}
