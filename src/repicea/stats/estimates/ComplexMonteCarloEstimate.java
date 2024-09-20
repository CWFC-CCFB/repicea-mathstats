/*
 * This file is part of the repicea-mathstats library.
 *
 * Copyright (C) 2023 His Majesty the King in Right of Canada
 * Author: Mathieu Fortin, Canadian Forest Service
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
package repicea.stats.estimates;

import java.security.InvalidParameterException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

import repicea.math.ComplexMatrix;
import repicea.math.ComplexSymmetricMatrix;
import repicea.math.HermitianMatrix;
import repicea.math.Matrix;
import repicea.stats.distributions.ComplexEmpiricalDistribution;

/**
 * An implementation of the Monte Carlo estimator for complex numbers.
 * @author Mathieu Fortin - November 2023
 */
@SuppressWarnings("serial")
public class ComplexMonteCarloEstimate extends ResamplingBasedEstimate<ComplexMatrix, HermitianMatrix>{

	/**
	 * Constructor.
	 */
	public ComplexMonteCarloEstimate() {
		super(new ComplexEmpiricalDistribution());
	}

	@Override
	public ComplexEmpiricalDistribution getDistribution() {
		return (ComplexEmpiricalDistribution) super.getDistribution();
	}
	
	/**
	 * Produce the real part of variance of diagonal element of the Hermitian matrix.
	 * @see ComplexEmpiricalDistribution#getVarianceRealPart()
	 * @return a Matrix instance
	 */
	public Matrix getVarianceRealPart() {
		return getDistribution().getVarianceRealPart();
	}

	/**
	 * Produce the imaginary part of variance of diagonal element of the Hermitian matrix.
	 * @see ComplexEmpiricalDistribution#getVarianceImaginaryPart()
	 * @return a Matrix instance
	 */
	public Matrix getVarianceImaginaryPart() {
		return getDistribution().getVarianceImaginaryPart();
	}

	/**
	 * Produce the pseudo-variance of the vector of complex random variables.
	 * @return a ComplexSymmetricMatrix instance
	 */
	public ComplexSymmetricMatrix getPseudoVariance() {
		return getDistribution().getPseudoVariance();
	}

//	/**
//	 * To be documented
//	 * @return
//	 */
//	public ComplexMonteCarloEstimate getInversedRealizations() {
//		ComplexMonteCarloEstimate est = new ComplexMonteCarloEstimate();
//		for (ComplexMatrix m : getDistribution().getInversedRealizations()) {
//			est.addRealization(m);
//		}
//		return est;
//	}

//	private class ComplexMonteCarloRealization implements Comparable<ComplexMonteCarloRealization> {
//		final double real;
//		final double absValue;
//		
//		ComplexMonteCarloRealization(ComplexNumber cn) {
//			real = cn.realPart;
//			absValue = cn.getAbsoluteValue();
//		}
//
//
//		@Override
//		public int compareTo(ComplexMonteCarloRealization o) {
//			if (absValue < o.absValue) {
//				return -1;
//			} else if (absValue > o.absValue) {
//				return 1;
//			} else {
//				return 0;
//			}
//		}
//		
//	}
	
	/**
	 * {@inheritDoc}<p>
	 * For this class, the realizations are ranked according to the value of their real part.
	 */
	@Override
	protected Matrix getQuantileForProbability(double probability) {
		if (probability < 0 || probability > 1) {
			throw new InvalidParameterException("The percentile must be between 0 and 1!");
		}
		List<ComplexMatrix> realizations = getRealizations();
		List<Double> realizationsForThisRow;
		int nbRows = realizations.get(0).m_iRows;
		Matrix percentileValues = new Matrix(nbRows,1);
		for (int i = 0; i < nbRows; i++) {
			realizationsForThisRow = new ArrayList<Double>();
			for (int j = 0; j < realizations.size(); j++) { 
				realizationsForThisRow.add(realizations.get(j).getValueAt(i, 0).realPart);
			}
			Collections.sort(realizationsForThisRow);
			int index = (int) Math.round(probability * realizations.size()) - 1;
			if (index < 0) {
				index = 0;
			} 
			percentileValues.setValueAt(i, 0, realizationsForThisRow.get(index));
		}
//		List<ComplexMonteCarloRealization> realizationsForThisRow;
//		int nbRows = realizations.get(0).m_iRows;
//		Matrix percentileValues = new Matrix(nbRows,1);
//		for (int i = 0; i < nbRows; i++) {
//			realizationsForThisRow = new ArrayList<ComplexMonteCarloRealization>();
//			for (int j = 0; j < realizations.size(); j++) { 
//				realizationsForThisRow.add(new ComplexMonteCarloRealization(realizations.get(j).getValueAt(i, 0)));
//			}
//			Collections.sort(realizationsForThisRow);
//			int index = (int) Math.round(probability * realizations.size()) - 1;
//			if (index < 0) {
//				index = 0;
//			} 
//			percentileValues.setValueAt(i, 0, realizationsForThisRow.get(index).real);
//		}
		return percentileValues;
	}
	
}
