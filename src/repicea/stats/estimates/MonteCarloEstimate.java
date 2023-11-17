/*
 * This file is part of the repicea library.
 *
 * Copyright (C) 2009-2012 Mathieu Fortin for Rouge-Epicea
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

import repicea.math.Matrix;
import repicea.stats.distributions.EmpiricalDistribution;
import repicea.util.REpiceaTranslator;
import repicea.util.REpiceaTranslator.TextableEnum;

/**
 * This estimate contains the realizations of a Monte Carlo simulations.
 * @author Mathieu Fortin - October 2011
 */
public class MonteCarloEstimate extends ResamplingBasedEstimate<Matrix> {
	
	private static final long serialVersionUID = 20110912L;
	
	public enum MessageID implements TextableEnum {
		Mean("Mean", "Moyenne"),
		Lower("Lower", "Inf"),
		Upper("Upper", "Sup"),
		ProbabilityLevel("Probability level", "Niveau de probabilit\u00E9");

		MessageID(String englishText, String frenchText) {
			setText(englishText, frenchText);
		}
		
		@Override
		public void setText(String englishText, String frenchText) {
			REpiceaTranslator.setString(this, englishText, frenchText);
		}
		
		@Override
		public String toString() {return REpiceaTranslator.getString(this);}
		
	}

	/**
	 * Constructor.
	 */
	public MonteCarloEstimate() {
		super(new EmpiricalDistribution());
	}
	
	/**
	 * This method returns a MonteCarloEstimate instance that results from the subtraction of two 
	 * MonteCarloEstimate instances with the same number of realizations. 
	 * @param estimate2 the estimate that is subtracted from this estimate
	 * @return a MonteCarloEstimate instance
	 */
	protected MonteCarloEstimate subtract(MonteCarloEstimate estimate2) {
		if (getNumberOfRealizations() != estimate2.getNumberOfRealizations()) {
			throw new InvalidParameterException("The number of realizations is not consistent!");
		}
		MonteCarloEstimate outputEstimate = new MonteCarloEstimate();
		for (int i = 0; i < getNumberOfRealizations(); i++) {
			outputEstimate.addRealization(getRealizations().get(i).subtract(estimate2.getRealizations().get(i)));
		}
		return outputEstimate;
	}
	
	/**
	 * This method returns a MonteCarloEstimate instance that results from the sum of two 
	 * MonteCarloEstimate instances with the same number of realizations. 
	 * @param estimate2 the estimate that is added to this estimate
	 * @return a MonteCarloEstimate instance
	 */
	protected MonteCarloEstimate add(MonteCarloEstimate estimate2) {
		if (getNumberOfRealizations() != estimate2.getNumberOfRealizations()) {
			throw new InvalidParameterException("The number of realizations is not consistent!");
		}
		MonteCarloEstimate outputEstimate = new MonteCarloEstimate();
		for (int i = 0; i < getNumberOfRealizations(); i++) {
			outputEstimate.addRealization(getRealizations().get(i).add(estimate2.getRealizations().get(i)));
		}
		return outputEstimate;
	}

	/**
	 * This method returns a MonteCarloEstimate instance that results from the product of original 
	 * MonteCarloEstimate instance and a scalar. 
	 * @param scalar the multiplication factor
	 * @return a MonteCarloEstimate instance
	 */
	protected MonteCarloEstimate multiply(double scalar) {
		MonteCarloEstimate outputEstimate = new MonteCarloEstimate();
		for (int i = 0; i < getNumberOfRealizations(); i++) {
			outputEstimate.addRealization(getRealizations().get(i).scalarMultiply(scalar));
		}
		return outputEstimate;
	}

	
	/**
	 * This method returns the percentile of the Monte Carlo simulated distribution.
	 * @param percentile a value between 0 and 1
	 * @return a Matrix instance that contains the percentile values
	 */
	protected Matrix getQuantileForProbability(double percentile) {
		if (percentile < 0 || percentile > 1) {
			throw new InvalidParameterException("The percentile must be between 0 and 1!");
		}
		List<Matrix> realizations = getRealizations();
		List<Double> realizationsForThisRow;
		int nbRows = realizations.get(0).m_iRows;
		Matrix percentileValues = new Matrix(nbRows,1);
		for (int i = 0; i < nbRows; i++) {
			realizationsForThisRow = new ArrayList<Double>();
			for (int j = 0; j < realizations.size(); j++) {
				realizationsForThisRow.add(realizations.get(j).getValueAt(i, 0));
			}
			Collections.sort(realizationsForThisRow);
			int index = (int) Math.round(percentile * realizations.size()) - 1;
			if (index < 0) {
				index = 0;
			} 
			percentileValues.setValueAt(i, 0, realizationsForThisRow.get(index));
		}
		return percentileValues;
	}
	

	/**
	 * This method returns a subset of the MonteCarloEstimate. For instance, if the estimate is multivariate, it is then possible
	 * to extract a MonteCarloEstimate with only the first and second variates. 
	 * @param indices a List of Integer that are the indices of the variates to be extracted.
	 * @return a MonteCarloEstimate
	 */
	public MonteCarloEstimate extractSubEstimate(List<Integer> indices) {
		MonteCarloEstimate subEstimate = new MonteCarloEstimate();
		for (Matrix realization : getRealizations()) {
			if (realization.isColumnVector()) {
				subEstimate.addRealization(realization.getSubMatrix(indices, null));
			} else {
				subEstimate.addRealization(realization.getSubMatrix(null, indices));
			}
		}
		return subEstimate;
	}
	
	@Override
	protected boolean isMergeableEstimate(Estimate<?,?> estimate) {
		if (estimate instanceof MonteCarloEstimate) {
			if (((MonteCarloEstimate) estimate).getNumberOfRealizations() == getNumberOfRealizations()) {
				return true;
			};
		}
		return false;
	}
	
	@Override
	public String toString() {
		return "Monte Carlo estimate (mean = " + getMean() + ", n = " + getNumberOfRealizations();
	}

	@Override
	public ConfidenceInterval getConfidenceIntervalBounds(double oneMinusAlpha) {
		Matrix lowerBoundValue = getQuantileForProbability(.5 * (1d - oneMinusAlpha));
		Matrix upperBoundValue = getQuantileForProbability(1d - .5 * (1d - oneMinusAlpha));
		return new ConfidenceInterval(lowerBoundValue, upperBoundValue, oneMinusAlpha);
	}
	
	
	@Override
	public Estimate<?,?> getDifferenceEstimate(Estimate<?,?> estimate2) {
		if (this.isMergeableEstimate(estimate2)) {
			return subtract((MonteCarloEstimate) estimate2);
		} else {
			return super.getDifferenceEstimate(estimate2);
		}
	}

	@Override
	public Estimate<?,?> getSumEstimate(Estimate<?,?> estimate2) {
		if (this.isMergeableEstimate(estimate2)) {
			return add((MonteCarloEstimate) estimate2);
		} else {
			return super.getSumEstimate(estimate2);
		}
	}

	@Override
	public Estimate<?,?> getProductEstimate(double scalar) {
		return multiply(scalar);
	}

	
	
}
