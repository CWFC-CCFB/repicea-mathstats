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
package repicea.stats.model.glm;

import repicea.math.Matrix;
import repicea.math.SymmetricMatrix;
import repicea.math.utility.GammaUtility;
import repicea.stats.distributions.utility.NegativeBinomialUtility;
import repicea.stats.model.glm.GeneralizedLinearModel.GLMIndividualLikelihood;

/**
 * This class simply handles the exponent y and 1-y of the likelihood function.  
 * @author Mathieu Fortin - October 2022
 */
@SuppressWarnings("serial")
public final class NegativeBinomialIndividualLikelihood extends GLMIndividualLikelihood {

	protected final LinkFunction linkFunction;
	private Matrix scale = new Matrix(1,1);
	
	public NegativeBinomialIndividualLikelihood(LinkFunction linkFunction) {
		super(linkFunction.getOriginalFunction());
		this.linkFunction = linkFunction;
	}

	@Override
	public void setParameters(Matrix beta) {
		for (int i = 0; i < beta.m_iRows; i++) {
			setParameterValue(i, beta.getValueAt(i, 0));
		}
	}

	@Override
	public Matrix getParameters() {
		Matrix beta = new Matrix(getNumberOfParameters(), 1);
		for (int i = 0; i < getNumberOfParameters(); i++) {
			beta.setValueAt(i, 0, getParameterValue(i));
		}
		return beta;
	}
	
	@Override
	public void setParameterValue(int index, double value) {
		if (additionalParameterIndices.contains(index)) {
			scale.setValueAt(0, 0, value);
		} else {
			super.setParameterValue(index, value);
		}
	}
	
	
	@Override
	public int getNumberOfParameters() {
		return super.getNumberOfParameters() + 1;
	}
	
	@Override
	public Matrix getPredictionVector() {
		Matrix mat = new Matrix(1,1);
		mat.setValueAt(0, 0, getMu());
		return mat;
	}
		
	@Override
	public Double getValue() {
		return NegativeBinomialUtility.getMassProbability(((Double) observedValues.getValueAt(0, 0)).intValue(), 
				getMu(), 
				getDispersion());
	}

	private double getMu() {
		return linkFunction.getValue();
	}
	
	private double getDispersion() {
		return scale.getValueAt(0,0);
	}
	
	@Override
	public double getParameterValue(int index) {
		if (additionalParameterIndices.contains(index)) {
			return scale.getValueAt(0, 0);
		} else {
			return super.getParameterValue(index);
		}
	}
	
	@Override
	public Matrix getGradient() {
		Matrix lfGradient = linkFunction.getGradient();
		double mu = getMu();
		double disp = getDispersion();
		double mu_x_disp = mu * disp;
		double invDisp = 1d / disp;
		int y = ((Double) observedValues.getValueAt(0, 0)).intValue();
		double f = getValue();
		double df_dmu =  f *  (y- mu)/(mu_x_disp * mu + mu);
		Matrix df_ddisp = new Matrix(1,1);
		df_ddisp.setValueAt(0, 0, 
				f * (
						(GammaUtility.digamma(invDisp) - GammaUtility.digamma(y + invDisp) + Math.log(mu_x_disp + 1))/ (disp * disp) + 
						(y - mu)/((mu_x_disp + 1)* disp))
					);
		return lfGradient.scalarMultiply(df_dmu).matrixStack(df_ddisp, true);
	}

	@Override
	public SymmetricMatrix getHessian() {
		SymmetricMatrix lfHessian = linkFunction.getHessian();
		Matrix lfGradient = linkFunction.getGradient();
		Matrix gradientPart = lfGradient.multiply(lfGradient.transpose());
		double f = getValue();
		double mu = getMu();
		double theta = getDispersion();
		double muTheta = mu * theta;
		double invTheta = 1d / theta;
		double y = observedValues.getValueAt(0, 0);
		double muThetaPlusOne = muTheta + 1;
		double digamma_invTheta = GammaUtility.digamma(invTheta);
		double digamma_yPlusInvTheta = GammaUtility.digamma(y + invTheta);
		
		double expTmp = (y/mu - (theta*y + 1) / muThetaPlusOne);
		double expTmp2 = (-y/(mu*mu) + theta*(theta*y + 1)/(muThetaPlusOne*muThetaPlusOne));
		gradientPart = gradientPart.scalarMultiply(f * (expTmp * expTmp + expTmp2));
		Matrix hessianPart = lfHessian.scalarMultiply(f * expTmp);
		Matrix d2f_d2mu = hessianPart.add(gradientPart);
		double df_dmu_part = (y - mu)/(muTheta*mu + mu);
		double df_dtheta_part = (digamma_invTheta - digamma_yPlusInvTheta + Math.log(muTheta + 1)) / 
				(theta*theta) + (y - mu)/((muTheta + 1)*theta);
		Matrix d2f_dtheta_dmu = lfGradient.scalarMultiply(f*(df_dmu_part*df_dtheta_part - (y - mu)/(muThetaPlusOne*muThetaPlusOne)));
				
		double y2 = y*y;
		double theta2 = theta*theta;
		double theta3 = theta2*theta;
		double mu2 = mu*mu;
		double logMuThetaPlusOne = Math.log(muThetaPlusOne);
		double trigamma_invTheta = GammaUtility.trigamma(1d/theta);
		double trigamma_yPlusInvTheta = GammaUtility.trigamma(y + 1d/theta);
		
		double longPoly = f /(theta2*theta2*muThetaPlusOne*muThetaPlusOne) * (
				3*theta3*mu2 - 2*theta3*mu2*logMuThetaPlusOne + theta2*mu2 + theta2*mu2*logMuThetaPlusOne*logMuThetaPlusOne
				- 2*theta2*mu2*logMuThetaPlusOne - theta2*mu2*trigamma_invTheta + 2*theta2*mu - 4*theta2*mu*logMuThetaPlusOne
				+ 2*theta*mu*logMuThetaPlusOne*logMuThetaPlusOne + logMuThetaPlusOne*logMuThetaPlusOne - 2*theta*mu*logMuThetaPlusOne - 2*theta*logMuThetaPlusOne
				- 2*theta*mu*trigamma_invTheta + muThetaPlusOne*muThetaPlusOne*digamma_invTheta*digamma_invTheta - trigamma_invTheta + theta2*y2 - 2*theta3*mu*y 
				+ theta2*mu2*trigamma_yPlusInvTheta - 2*theta2*mu*y + 2*theta2*mu*y*logMuThetaPlusOne - theta2*y
				+ 2*theta*y*logMuThetaPlusOne + 2*theta*mu*trigamma_yPlusInvTheta + muThetaPlusOne*muThetaPlusOne*digamma_yPlusInvTheta*digamma_yPlusInvTheta
				- 2*muThetaPlusOne*digamma_invTheta*(theta*(theta*mu + mu - y +1) - muThetaPlusOne*logMuThetaPlusOne)
				+ 2*muThetaPlusOne*digamma_yPlusInvTheta*(-muThetaPlusOne*logMuThetaPlusOne - muThetaPlusOne*digamma_invTheta + theta*(theta*mu + mu - y + 1))
				+ trigamma_yPlusInvTheta);
		
		Matrix d2f_d2theta = new Matrix(1,1,longPoly,0);
				
		Matrix hessian = d2f_d2mu.matrixStack(d2f_dtheta_dmu.transpose(), true)
				.matrixStack(d2f_dtheta_dmu.matrixStack(d2f_d2theta, true), false);
		return SymmetricMatrix.convertToSymmetricIfPossible(hessian);	
	}


}
