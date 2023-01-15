/*
 * This file is part of the repicea library.
 *
 * Copyright (C) 2009-2022 Mathieu Fortin for Rouge-Epicea
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
package repicea.stats.model.dist;

import java.security.InvalidParameterException;
import java.util.ArrayList;
import java.util.List;

import repicea.math.LogFunctionWrapper;
import repicea.math.Matrix;
import repicea.math.functions.GaussianFunction;
import repicea.stats.estimators.Estimator;
import repicea.stats.estimators.MaximumLikelihoodEstimator;
import repicea.stats.estimators.MaximumLikelihoodEstimator.MaximumLikelihoodCompatibleModel;
import repicea.stats.model.AbstractStatisticalModel;
import repicea.stats.model.CompositeLogLikelihood;
import repicea.stats.model.IndividualLogLikelihood;
import repicea.stats.model.SimpleCompositeLogLikelihood;

/**
 * The GaussianModel class makes it possible to fit a Gaussian distribution to a list of values.<br>
 * <br>
 * The fit relies on a maximum likelihood estimator.
 * 
 * @author Mathieu Fortin - July 2022
 */
public class GaussianModel extends AbstractStatisticalModel implements MaximumLikelihoodCompatibleModel {

	@SuppressWarnings("serial")
	private class GaussianLogLikehood extends LogFunctionWrapper implements IndividualLogLikelihood {

		private GaussianLogLikehood() {
			super(new GaussianFunction());
		}

		@Override
		public void setYVector(Matrix yVector) {
			if (yVector.getNumberOfElements() != 1) {
				throw new InvalidParameterException("The yVector should be a unique element!");
			}
			getOriginalFunction().setVariableValue(0, yVector.getValueAt(0, 0));
		}

		@Override
		public Matrix getYVector() {
			return new Matrix(1, 1, getOriginalFunction().getVariableValue(0), 0d);
		}

		@Override
		public Matrix getPredictionVector() {return null;}
	}

	
	private final List<Double> values;
	private final SimpleCompositeLogLikelihood cLL;
	private final IndividualLogLikelihood individualLLK;

	/**
	 * General constructor.
	 * @param values a sample of the distribution
	 * @param startingValues a 2x1 matrix with starting values. If set to null, the starting values are then 0 and 1 
	 * for mu and sigma2, respectively.
	 */
	public GaussianModel(List<Double> values, Matrix startingValues) {
		super();
		this.values = new ArrayList<Double>();
		this.values.addAll(values);
		this.individualLLK = new GaussianLogLikehood();
		cLL = new SimpleCompositeLogLikelihood(individualLLK, new Matrix(values));
		setParameters(startingValues);
		try {
			setModelDefinition("pdf(y) = 1/(2*PI*sigma2)^(1/2) * e^(-(y-mu)^2 / (2 * sigma2))");
		} catch (Exception e) {}
	}
	
	/**
	 * Constructor based on default starting values. <br>
	 * <br>
	 * The mu and sigma2 parameters are set to 0 and 1, respectively.
	 * @param values a sample of the distribution
	 */
	public GaussianModel(List<Double> values) {
		this(values, null);
	}
	
	@Override
	public void setParameters(Matrix beta) {
		if (beta == null) {
			Matrix betaDefault = new Matrix(2,1);
			betaDefault.setValueAt(1, 0, 1d);
			individualLLK.setParameters(betaDefault);
		} else {
			individualLLK.setParameters(beta);
		}
	}

//	@Override
//	public Matrix getParameters() {
//		return individualLLK.getParameters();
//	}

	protected Estimator instantiateDefaultEstimator() {return new MaximumLikelihoodEstimator(this);}

	@Override
	public boolean isInterceptModel() {return false;}

	@Override
	public List<String> getEffectList() {
		List<String> effectList = new ArrayList<String>(); 
		effectList.add("mu parameter");
		effectList.add("sigma2 parameter");
		return effectList;
	}

	@Override
	public int getNumberOfObservations() {return values.size();}

	@Override
	public double getConvergenceCriterion() {
		return 1E-8;
	}

	@Override
	public CompositeLogLikelihood getCompleteLogLikelihood() {return cLL;}
	
	@Override
	public String toString() {
		return "Gaussian model";
	}

	@Override
	public List<String> getOtherParameterNames() {return new ArrayList<String>();}

}
