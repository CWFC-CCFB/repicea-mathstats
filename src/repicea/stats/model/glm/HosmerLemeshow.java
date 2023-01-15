package repicea.stats.model.glm;

import java.security.InvalidParameterException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

import repicea.math.Matrix;

public class HosmerLemeshow {

	@SuppressWarnings("rawtypes")
	protected class Couplet implements Comparable {

		private double yValue;
		private double predicted;
		
		protected Couplet(double yValue, double predicted) {
			this.yValue = yValue;
			this.predicted = predicted;
		}
		
		
		@Override
		public int compareTo(Object o) {
			Couplet couplet = (Couplet) o;
			if (predicted < couplet.predicted) {
				return -1;
			} else if (predicted > couplet.predicted) {
				return 1;
			} else {
				return 0;
			}
		}
	}
	
	
	@SuppressWarnings("serial")
	protected final class HosmerLemeshowSubgroup extends ArrayList<Couplet> {
		
		protected HosmerLemeshowSubgroup() {}
		
		private double getTotalPositiveResponse() {
			double sum = 0;
			for (int i = 0; i < size(); i++) {
				sum += get(i).yValue;
			}
			return sum;
		}

		private double getAveragePredicted() {
			double sum = 0;
			for (int i = 0; i < size(); i++) {
				sum += get(i).predicted;
			}
			return sum / size();
		}

		
		protected double getValue() {
			double averagePred = getAveragePredicted();
			double totalPred = size() * averagePred;
			double diff = (getTotalPositiveResponse() - totalPred);
			double value = diff * diff / (totalPred * (1 - averagePred));
			return value;
		}
				
		
	}

	
	private final GeneralizedLinearModel glm;
	private final List<HosmerLemeshowSubgroup> subgroups;
	private final int numberOfGroups;
	private double cValue;
	

	/**
	 * General constructor.
	 * @param glm
	 * @param numberOfGroups
	 */
	public HosmerLemeshow(GeneralizedLinearModel glm, int numberOfGroups) {
		if (numberOfGroups < 2) {
			throw new InvalidParameterException("The number of groups must be larger than 2!");
		}
		scanResponseVariable(glm);
		this.numberOfGroups = numberOfGroups;
		this.glm = glm;
		subgroups = new ArrayList<HosmerLemeshowSubgroup>();
		computeCValue();
	}

	
	
	
	
	
	
	private void scanResponseVariable(GeneralizedLinearModel glm) {
		Matrix responseVector = glm.y;
		for (int i = 0; i < responseVector.m_iRows; i++) {
			if (responseVector.getValueAt(i, 0) != 0d && responseVector.getValueAt(i, 0) != 1d) {
				throw new InvalidParameterException("The response variable must be either 0 or 1!");
			}
		}
	}


	/**
	 * Constructor with number of groups set to 10.
	 * @param glm a GeneralizedLinearModel instance
	 */
	public HosmerLemeshow(GeneralizedLinearModel glm) {
		this(glm, 10);
	}
	

	@SuppressWarnings("unchecked")
	private void computeCValue() {
		List<Couplet> couplets = new ArrayList<Couplet>();
		Matrix observed = glm.y;
		Matrix predicted = glm.getPredicted();
		for (int i = 0; i < observed.m_iRows; i++) {
			couplets.add(new Couplet(observed.getValueAt(i, 0), predicted.getValueAt(i, 0)));
		}
		
		Collections.sort(couplets);
		
		int numberObservationsByGroup = (int) Math.floor(couplets.size() * 1d / numberOfGroups);
		int diff = couplets.size() - numberObservationsByGroup * 10; 
		List<Integer> groupSizes = new ArrayList<Integer>(); 
		for (int i = 0; i < numberOfGroups; i++) {
			if (diff > 0) {
				diff--;
				groupSizes.add(numberObservationsByGroup + 1);
			} else {
				groupSizes.add(numberObservationsByGroup);
			}
		}

		HosmerLemeshowSubgroup subgroup;
		int index = 0;
		for (int g = 0; g < groupSizes.size(); g++) {
			subgroup = new HosmerLemeshowSubgroup();
			subgroups.add(subgroup);
			for (int i = index; i < index + groupSizes.get(g); i++) {
				subgroup.add(couplets.get(i));
			}
			index += groupSizes.get(g);
		}
		
		double stat = 0;
		for (int g = 0; g < subgroups.size(); g++) {
			stat += subgroups.get(g).getValue();
		}
		
		cValue = stat;
	}
	
	@Override
	public String toString() {
		return "C Statistic = " + cValue + ", with " + (numberOfGroups - 2) + " degrees of freedom"; 
	}
	
}
