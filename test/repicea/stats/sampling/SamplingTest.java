package repicea.stats.sampling;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.junit.Assert;
import org.junit.Test;


public class SamplingTest {

	@SuppressWarnings({ "rawtypes", "unchecked" })
	@Test
	public void pickARandomSamplingWithoutReplacement() {
		List population = new ArrayList();
		for (int i = 0; i < 100; i++) {
			population.add(i);
		}
		
		List sample = SamplingUtility.getSample(population, 50);
		
		Map frequencyMap = SamplingUtility.getObservationFrequencies(sample);
		
		Assert.assertEquals("Testing if it was without replacement", sample.size(), frequencyMap.size()); 
	}
	
	@SuppressWarnings({ "rawtypes", "unchecked" })
	@Test
	public void pickARandomSamplingWithReplacement() {
		List population = new ArrayList();
		for (int i = 0; i < 100; i++) {
			population.add(i);
		}
		
		List sample = SamplingUtility.getSample(population, 100, true);
		
		Map frequencyMap = SamplingUtility.getObservationFrequencies(sample);
		
		Assert.assertEquals("Testing if it was with replacement", sample.size() > frequencyMap.size(), true); 
	}

	
	@SuppressWarnings({ "rawtypes", "unchecked" })
	@Test
	public void pickARandomSamplingWithoutReplacementFromMap() {
		Map population = new HashMap();
		for (int i = 0; i < 100; i++) {
			population.put(i,i);
		}
		
		Map sample = SamplingUtility.getSample(population, 50);
		
		Map frequencyMap = SamplingUtility.getObservationFrequencies(sample.values());
		
		Assert.assertEquals("Testing if it was without replacement", sample.size(), frequencyMap.size()); 
	}
	

	
}
