package repicea.stats.estimates;

import java.security.InvalidParameterException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.concurrent.ConcurrentHashMap;

import repicea.math.Matrix;
import repicea.stats.sampling.PopulationUnit;

@SuppressWarnings("serial")
public abstract class AbstractSimplePointEstimate extends AbstractPointEstimate {

	private final Map<String, PopulationUnit> observations;
	private final double populationSize;

	
	/**
	 * Basic constructor without population size.
	 */
	protected AbstractSimplePointEstimate() {
		super();
		observations = new ConcurrentHashMap<String, PopulationUnit>();
		populationSize = -1d;
	}

	/**
	 * Constructor with population size.
	 * @param populationSize the number of units in the population.
	 */
	protected AbstractSimplePointEstimate(double populationSize) {
		super();
		if (populationSize <= 0) {
			throw new InvalidParameterException("The population size must be greater than 0!");
		}
		observations = new ConcurrentHashMap<String, PopulationUnit>();
		this.populationSize = populationSize;
	}
	

	/**
	 * Create a Matrix instance with each row representing one observation. The order is ensured by
	 * the list of sample Ids.
	 * @return a Matrix
	 */
	protected Matrix getObservationMatrix() {
		Matrix outputMatrix = null;
		int nbObservations = getObservations().size();
		int nbElementsPerObs = 0;
		List<String> sampleIds = getPopulationUnitIds();
		for (int i = 0; i < sampleIds.size(); i++) {
			String sampleId = sampleIds.get(i);
			PopulationUnit obs = getObservations().get(sampleId);
			if (outputMatrix == null) {
				nbElementsPerObs = obs.getData().m_iRows;
				outputMatrix = new Matrix(nbObservations, nbElementsPerObs);
			}
			outputMatrix.setSubMatrix(obs.getData().transpose(), i, 0);
		}
		return outputMatrix;
	}

	/**
	 * Create a List with the ordered sample ids 
	 * @return a List instance
	 */
	@Override
	protected final List<String> getPopulationUnitIds() {
		List<String> puIds = new ArrayList<String>();
		puIds.addAll(observations.keySet());
		Collections.sort(puIds);
		return puIds;
	}
	

	/**
	 * Add an observation to the sample.
	 * 
	 * @param obs a PopulationUnitObservation instance
	 */
	public void addObservation(PopulationUnit obs) {
		if (obs == null) {
			throw new InvalidParameterException("The obs argument must be non null!");
		}
		String sampleId = obs.getSampleId();
		if (observations.containsKey(sampleId)) {
			throw new InvalidParameterException("The sample id " + sampleId + " is already contained in the observation map!");
		}
		if (nCols == 0) {
			nCols = obs.getData().m_iCols;
		}
		if (nRows == 0) {
			nRows = obs.getData().m_iRows;
		}
		if (obs.getData().m_iCols != nCols || obs.getData().m_iRows != nRows) {
			throw new InvalidParameterException("The observation is incompatible with what was already observed!");
		} else {
			observations.put(sampleId, obs);
		}
	}

	protected final Map<String, PopulationUnit> getObservations() {return observations;}

	@Override
	public double getPopulationSize() {return populationSize;}
	

	@Override
	protected boolean isMergeableEstimate(Estimate<?,?,?> estimate) {
		if (getClass().equals(estimate.getClass())) {
			AbstractSimplePointEstimate pe = (AbstractSimplePointEstimate) estimate;
			if (getPopulationUnitIds().equals(pe.getPopulationUnitIds())) {	// make sure we have the same sample ids
				if (nRows == pe.nRows) {
					if (nCols == pe.nCols) {
						if (getPopulationSize() == pe.getPopulationSize()) {
							return true;
						}
					}
				}
			}
		}
		return false;
	}

	@Override
	protected abstract AbstractSimplePointEstimate add(PointEstimate pointEstimate);

	@Override
	protected abstract AbstractSimplePointEstimate subtract(PointEstimate pointEstimate);

	@Override
	protected abstract AbstractSimplePointEstimate multiply(double scalar);

	@Override
	public final int getSampleSize() {return observations.size();}
}
