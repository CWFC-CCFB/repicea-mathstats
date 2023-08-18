/*
 * This file is part of the repicea-mathstats library.
 *
 * Copyright (C) 2023 His Majesty the King in Right of Canada
 * Author: Mathieu Fortin - Canadian Wood Fibre Centre, CFS
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
package repicea.stats.model.lm;

import java.security.InvalidParameterException;

import repicea.math.AbstractMathematicalFunctionWrapper;
import repicea.math.Matrix;
import repicea.math.SymmetricMatrix;
import repicea.stats.LinearStatisticalExpression;
import repicea.stats.model.IndividualLogLikelihood;

@SuppressWarnings("serial")
public class GaussianIndividualLogLikelihood extends AbstractMathematicalFunctionWrapper implements IndividualLogLikelihood {

	protected Matrix yVector;
	protected double sigma2;
	protected final int nbParmsFromMatrixX;
	
	protected GaussianIndividualLogLikelihood(int nbParmsFromMatrixX, Matrix startingValues) {
		super(new LinearStatisticalExpression());
		this.nbParmsFromMatrixX = nbParmsFromMatrixX;
		getOriginalFunction().setParameters(startingValues.getSubMatrix(0, nbParmsFromMatrixX - 1, 0, 0));
		sigma2 = startingValues.getValueAt(startingValues.m_iRows - 1, 0);
	}
	
	@Override
	public void setYVector(Matrix yVector) {this.yVector = yVector;}

	@Override
	public Matrix getYVector() {return yVector;}

	@Override
	public Matrix getPredictionVector() {
		return new Matrix(1,1,getOriginalFunction().getValue(),0);
	}
	
	@Override
	public void setParameterValue(int i, double value) {
		if (i < nbParmsFromMatrixX) {
			super.setParameterValue(i, value);
		} else if (i == nbParmsFromMatrixX) {
			sigma2 = value;
		} else {
			throw new InvalidParameterException("There is no parameter defined by index " + i);
		}
	}
	
	@Override
	public void setParameters(Matrix beta) {
		for (int i = 0; i < beta.m_iRows; i++) {
			setParameterValue(i, beta.getValueAt(i, 0));
		}
	}
	
	@Override
	public double getParameterValue(int i) {
		if (i < nbParmsFromMatrixX) {
			return super.getParameterValue(i);
		} else if (i == nbParmsFromMatrixX) {
			return sigma2;
		} else {
			throw new InvalidParameterException("There is no parameter defined by index " + i);
		}
	}
	
	@Override
	public Matrix getParameters() {
		Matrix beta = new Matrix(getNumberOfParameters(), 1);
		for (int i = 0; i < beta.m_iRows; i++) {
			beta.setValueAt(i, 0, getParameterValue(i));
		}
		return beta;
	}

	@Override
	public int getNumberOfParameters() {
		return getOriginalFunction().getNumberOfParameters() + 1;
	}

	@Override
	public Double getValue() {
		double res = yVector.getValueAt(0, 0) - getOriginalFunction().getValue();
		return - 0.5 * (res * res) / sigma2 - 0.5 * Math.log(2 * Math.PI) - 0.5 * Math.log(sigma2);
	}

	
	@Override
	public Matrix getGradient() {
		Matrix g = new Matrix(getNumberOfParameters(), 1);
		double res = yVector.getValueAt(0, 0) - getOriginalFunction().getValue();
		double df_dMu = res / sigma2;
		Matrix dMu_dBeta = getOriginalFunction().getGradient(); 
		g.setSubMatrix(dMu_dBeta.scalarMultiply(df_dMu), 0, 0);
		double df_dsigma2 = 0.5 * (res * res) / (sigma2 * sigma2) - 0.5 / sigma2;
		g.setValueAt(g.m_iRows - 1, 0, df_dsigma2);
		return g;
	}

	@Override
	public SymmetricMatrix getHessian() {
		int nbParms = getNumberOfParameters();
		double res = yVector.getValueAt(0, 0) - getOriginalFunction().getValue();
		double df_dMu = res / sigma2;
		Matrix dMu_dBeta = getOriginalFunction().getGradient();
//		Matrix d2Mu_d2Beta = getOriginalFunction().getHessian();	// unnecessary since it is 0 by definition
//		Matrix hessianPart1_d2f_d2Beta = d2Mu_d2Beta.scalarMultiply(df_dMu);
		Matrix hessianPart2_d2f_d2Beta = dMu_dBeta.multiply(dMu_dBeta.transpose()).scalarMultiply(-1/sigma2);
		double df_dMu_dSigma2 = - df_dMu / sigma2;
		Matrix hessian_d2f_dBeta_dSigma2 = dMu_dBeta.scalarMultiply(df_dMu_dSigma2);
		double d2f_d2Sigma2 = - (res * res) / (sigma2 * sigma2 * sigma2) + 0.5 / (sigma2 * sigma2);
		Matrix hessian = new Matrix(nbParms, nbParms);
//		hessian.setSubMatrix(hessianPart1_d2f_d2Beta.add(hessianPart2_d2f_d2Beta), 0, 0);
		hessian.setSubMatrix(hessianPart2_d2f_d2Beta, 0, 0);
		hessian.setSubMatrix(hessian_d2f_dBeta_dSigma2, 0, hessian.m_iCols - 1);
		hessian.setSubMatrix(hessian_d2f_dBeta_dSigma2.transpose(), hessian.m_iRows - 1, 0);
		hessian.setValueAt(hessian.m_iRows - 1, hessian.m_iCols - 1, d2f_d2Sigma2);
		return SymmetricMatrix.convertToSymmetricIfPossible(hessian);
	}

	
	
}
