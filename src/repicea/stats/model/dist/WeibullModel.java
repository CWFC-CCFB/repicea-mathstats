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
import java.util.Collections;
import java.util.List;

import repicea.math.LogFunctionWrapper;
import repicea.math.Matrix;
import repicea.math.functions.WeibullFunction;
import repicea.stats.estimators.Estimator;
import repicea.stats.estimators.MaximumLikelihoodEstimator;
import repicea.stats.estimators.MaximumLikelihoodEstimator.MaximumLikelihoodCompatibleModel;
import repicea.stats.model.AbstractStatisticalModel;
import repicea.stats.model.CompositeLogLikelihood;
import repicea.stats.model.IndividualLogLikelihood;
import repicea.stats.model.SimpleCompositeLogLikelihood;

/**
 * Fit the Weibull distribution to a sample of an unknown distribution. <br>
 * <br>
 * The model is fitted using a maximum likelihood estimator and the Weibull distribution
 * is parameterized as follows: <br>
 * <br>
 * pdf(y) = k/lambda * ((y - theta)/lambda)^(k-1) * e^(-((y - theta)/lambda)^k) <br>
 * <br>
 * 
 * When not supplied the starting parameters are estimated using the percentile method shown in 
 * <a href=https://doi.org/10.1093/forestscience/31.1.260> Zarnoch, S. J. and Dell, T. R. 1985. An
 * evaluation of percentile and maximum likelihood estimators of Weibull parameters. Forest Science
 * 31(1): 260-268.</a> 
 * 
 * @author Mathieu Fortin - July 2022
 */
public class WeibullModel extends AbstractStatisticalModel implements MaximumLikelihoodCompatibleModel {

	
	@SuppressWarnings("serial")
	private class WeibullLogLikehood extends LogFunctionWrapper implements IndividualLogLikelihood {

		private WeibullLogLikehood(boolean isLocationParameterEnabled) {
			super(isLocationParameterEnabled ? new WeibullFunction(1,1,0) : new WeibullFunction(1,1));
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

	
	
//	private class WeibullLogLikehood extends AbstractMathematicalFunction implements IndividualLogLikelihood {
//
////		private final Matrix parameters;
//		private Matrix yVector;
//		
//		private WeibullLogLikehood() {
//			if (isLocationParameterEnabled) {
////				parameters = new Matrix(3, 1); 
//				setParameters(new Matrix(3, 1));
//			} else {
//				parameters = new Matrix(2, 1);
//			}
//		 }
//		
//		@Override
//		public void setYVector(Matrix yVector) {
//			if (yVector.m_iCols > 1 || yVector.m_iRows > 1)
//				throw new InvalidParameterException("The y value should be embedded in a 1x1 matrix!");
//			this.yVector = yVector;
//		}
//
//		@Override
//		public Matrix getYVector() {return yVector;}
//
//		@Override
//		public Matrix getPredictionVector() {
//			// TODO Auto-generated method stub
//			return null;
//		}
//
//		@Override
//		public int getNumberOfParameters() {return parameters.m_iRows;}
//
//		@Override
//		public int getNumberOfVariables() {return 0;}
//
//		@Override
//		public Matrix getParameters() {return parameters;}
//
//		@Override
//		public void setBounds(int parameterIndex, ParameterBound bound) {}
//
//		@Override
//		public void setVariableValue(int variableIndex, double variableValue) {
//			throw new InvalidParameterException("The WeibullModel class does not implement variables!");
//		}
//
//		@Override
//		public void setParameterValue(int parameterIndex, double parameterValue) {
//			parameters.setValueAt(parameterIndex, 0, parameterValue);
//		}
//
//		@Override
//		public void setParameters(Matrix beta) {
//			for (int i = 0; i < beta.m_iRows; i++) {
//				setParameterValue(i, beta.getValueAt(i, 0));
//			}
//		}
//
//		@Override
//		public void setVariables(Matrix xVector) {
//			throw new InvalidParameterException("The WeibullModel class does not implement variables!");
//		}
//
//		@Override
//		public double getVariableValue(int variableIndex) {
//			throw new InvalidParameterException("The WeibullModel class does not implement variables!");
//		}
//
//		@Override
//		public double getParameterValue(int parameterIndex) {
//			return parameters.getValueAt(parameterIndex, 0);
//		}
//		
//		@Override
//		public Double getValue() {
//			double k = parameters.getValueAt(0, 0);
//			double lambda = parameters.getValueAt(1, 0);
//			double yMod = isLocationParameterEnabled ? 
//					yVector.getValueAt(0, 0) - parameters.getValueAt(2, 0) :
//						yVector.getValueAt(0, 0);
//			double value = Math.log(k) + (k-1) * Math.log(yMod) - k * Math.log(lambda) - Math.pow(yMod / lambda, k);
//			if (Double.isNaN(value)) {
//				int u = 0;
//			}
//			return value;
//		}
//
//		@Override
//		public Matrix getGradient() {
//			double k = parameters.getValueAt(0, 0);
//			double lambda = parameters.getValueAt(1, 0);
//			double yMod = isLocationParameterEnabled ? 
//					yVector.getValueAt(0, 0) - parameters.getValueAt(2, 0) :
//						yVector.getValueAt(0, 0);
//			Matrix gradient = new Matrix(parameters.m_iRows, 1);
//			double yMod_lambda = yMod / lambda;
//			double yMod_lambda_powK = Math.pow(yMod_lambda, k);
//			double k_lambda = k/lambda;
//			
//			double dl_dk = 1/k + Math.log(yMod) - Math.log(lambda) - yMod_lambda_powK * Math.log(yMod_lambda);
// 			gradient.setValueAt(0, 0, dl_dk);
//   			
//			double dl_dlambda = -k_lambda + k_lambda * yMod_lambda_powK;
// 			gradient.setValueAt(1, 0, dl_dlambda);
//
// 			if (isLocationParameterEnabled) {
// 				double dl_dtheta = -(k-1)/yMod + k_lambda * Math.pow(yMod_lambda, k - 1); 
// 	 			gradient.setValueAt(2, 0, dl_dtheta);
// 			}
// 			
//			return gradient;
//		}
//
//		@Override
//		public Matrix getHessian() {
//			Matrix hessian = new Matrix(parameters.m_iRows, parameters.m_iRows);
//			double k = parameters.getValueAt(0, 0);
//			double lambda = parameters.getValueAt(1, 0);
//			double yMod = isLocationParameterEnabled ? 
//					yVector.getValueAt(0, 0) - parameters.getValueAt(2, 0) :
//						yVector.getValueAt(0, 0);
//			double yMod_lambda = yMod / lambda;
//			double yMod_lambda_powK = Math.pow(yMod_lambda, k);
//			double logYMod_lambda = Math.log(yMod_lambda);
//			double k_lambda = k/lambda;
//
//			double d2l_d2k = -1/(k*k) - yMod_lambda_powK * logYMod_lambda * logYMod_lambda;
//			hessian.setValueAt(0, 0, d2l_d2k);
//			
//			double d2l_d2lambda = k_lambda / lambda * (-yMod_lambda_powK + 1 - k *yMod_lambda_powK);
//			hessian.setValueAt(1, 1, d2l_d2lambda);
//
//			double d2l_dkdlambda = -1/lambda + yMod_lambda_powK / lambda * (k * logYMod_lambda + 1);
//			hessian.setValueAt(0, 1, d2l_dkdlambda);
//			hessian.setValueAt(1, 0, d2l_dkdlambda);
//			
//			if (isLocationParameterEnabled) {
//				double d2l_d2theta = -(k-1)/(yMod * yMod) - k_lambda/lambda * (k-1) * Math.pow(yMod_lambda, k-2); 
//				hessian.setValueAt(2, 2, d2l_d2theta);
//				
//				double d2l_dkdtheta = -1/yMod + Math.pow(yMod_lambda, k - 1) / lambda * (k * logYMod_lambda + 1);
//				hessian.setValueAt(2, 0, d2l_dkdtheta);
//				hessian.setValueAt(0, 2, d2l_dkdtheta);
//				
//				double d2l_dlambdadtheta = - k_lambda * k_lambda * Math.pow(yMod_lambda, k - 1);
//				hessian.setValueAt(2, 1, d2l_dlambdadtheta);
//				hessian.setValueAt(1, 2, d2l_dlambdadtheta);
//			}
//
//			return hessian;
//		}
//	}
		
	private final List<Double> values;
	private final SimpleCompositeLogLikelihood cLL;
	private final IndividualLogLikelihood individualLLK;
	private final boolean isLocationParameterEnabled;
	
	public WeibullModel(List<Double> values, boolean isLocationParameterEnabled, Matrix startingValues) {
		super();
		this.isLocationParameterEnabled = isLocationParameterEnabled;
		this.values = new ArrayList<Double>();
		this.values.addAll(values);
		this.individualLLK = new WeibullLogLikehood(isLocationParameterEnabled);
		cLL = new SimpleCompositeLogLikelihood(individualLLK, new Matrix(values));
		setParameters(startingValues);
		try {
			setModelDefinition(isLocationParameterEnabled ? 
					"pdf(y) = k/lambda * ((y - theta)/lambda)^(k-1) * e^(-((y - theta)/lambda)^k)" :
						"pdf(y) = k/lambda * (y/lambda)^(k-1) * e^(-(y/lambda)^k)");
		} catch (Exception e) {}
	}

	public WeibullModel(List<Double> values, boolean isLocationParameterEnabled) {
		this(values, isLocationParameterEnabled, null);
	}
	
	@Override
	public void setParameters(Matrix beta) {
		if (beta == null) {
			Matrix betaDefault = isLocationParameterEnabled ? new Matrix(3,1) : new Matrix(2,1);
			
			List<Double> valuesTmp = new ArrayList<Double>();
			valuesTmp.addAll(values);
			Collections.sort(valuesTmp);
			double theta_pct = 0d;
			if (isLocationParameterEnabled) {
				double x_1 = valuesTmp.get(0);
				double x_2 = valuesTmp.get(1);
				double x_n = valuesTmp.get(valuesTmp.size() - 1);
				theta_pct = (x_1 * x_n - x_2 * x_2) / (x_1 + x_n - 2 * x_2);
				if (theta_pct < 0) {
					theta_pct = 0;
				}
				betaDefault.setValueAt(2, 0, theta_pct);
			}
			
			double lambda_pct = -theta_pct + valuesTmp.get((int) (Math.round(0.63 * valuesTmp.size()) - 1));
			betaDefault.setValueAt(1, 0, lambda_pct);
			
			double p_k = 0.97366;
			double p_i = 0.16731;
			double numerator = Math.log(Math.log(1 - p_k) / Math.log(1 - p_i));
			int index_p_k = (int) Math.round(valuesTmp.size() * p_k) - 1;
			int index_p_i = (int) Math.round(valuesTmp.size() * p_i) - 1;
			double denominator = Math.log((valuesTmp.get(index_p_k) - theta_pct) / (valuesTmp.get(index_p_i) - theta_pct));
			double k_pct = numerator / denominator;
			betaDefault.setValueAt(0, 0, k_pct);
			
			individualLLK.setParameters(betaDefault);
		} else {
			individualLLK.setParameters(beta);
		}
	}


	@Override
	protected Estimator instantiateDefaultEstimator() {return new MaximumLikelihoodEstimator(this);}

	@Override
	public boolean isInterceptModel() {return false;}

	@Override
	public List<String> getEffectList() {
		List<String> effectList = new ArrayList<String>(); 
		effectList.add("Shape parameter");
		effectList.add("Scale parameter");
		if (isLocationParameterEnabled) {
			effectList.add("Location parameter");
		}
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
		return "Weibull model";
	}

	@Override
	public List<String> getOtherParameterNames() {return new ArrayList<String>();}

}
