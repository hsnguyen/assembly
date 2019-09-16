package org.rtassembly.npgraph;

import java.util.HashMap;
import java.util.Set;
import com.google.common.collect.Sets;
/**
 * Class to represent multiplicity of edge or node
 * in term of a list of bin:count key-value pairs
 * @author sonhoanghguyen
 *
 */
public class Multiplicity {
	HashMap<PopBin,Integer> mmap;
	
	public Multiplicity() {
		mmap=new HashMap<>();
	}
	public Multiplicity(HashMap<PopBin,Integer> map){
		this();
		for(PopBin b:map.keySet())
			if(map.get(b)>0)
				mmap.put(b, map.get(b));
	}
	
	public Multiplicity(PopBin b, int m){
		this();
		mmap.put(b,m);
	}
	
	public HashMap<PopBin, Integer> getFullMap(){
		return mmap;
	}
	
	public Set<PopBin> getBinsSet(){
		return mmap.keySet();
	}
	public int getBinCount(PopBin b){
		if(!mmap.containsKey(b))
			return 0;
		else
			return mmap.get(b);
	}
	public void add(Multiplicity m){
		for(PopBin b:m.getBinsSet())
			mmap.put(b,getBinCount(b) + m.getBinCount(b));
		
	}
	public double getCoverage(){
		double retval=0;
		for(PopBin b:getBinsSet())
			retval+=mmap.get(b)*b.estCov;
		return retval;
	}
	/**
	 * Substracting this to other Multiplicity
	 * this=minuend
	 * @param subtrahend
	 * @return new Multiplicity object = difference
	 */
	public Multiplicity substract(Multiplicity subtrahend){
		HashMap<PopBin, Integer> retval = new HashMap<PopBin, Integer>();
		for(PopBin b:Sets.union(getBinsSet(), subtrahend.getBinsSet())){
			int _minuend = getBinCount(b),
				_diff = _minuend - (subtrahend.getBinCount(b));	
			if(_diff >= 0)
				retval.put(b, _diff);
			else
				return null; //something wrong! TODO: make it right?
		}
		
		return new Multiplicity(retval);
	}
	/**
	 * Adding this and other Multiplicity
	 * this=augend
	 * @param addend
	 * @return new Multiplicity object = sum
	 */
	public Multiplicity sum(Multiplicity addend){
		HashMap<PopBin, Integer> retval = new HashMap<PopBin, Integer>();
		for(PopBin b:Sets.union(getBinsSet(), addend.getBinsSet())){
			int sum=getBinCount(b) + addend.getBinCount(b);
			if(sum>=0)
				retval.put(b, sum);
			else 
				return null; //something wrong
		}
		
		return new Multiplicity(retval);
	}
	
	//get sum of all bin counts
	public int getSum(){
		return mmap.values().stream().mapToInt(i->i.intValue()).sum();
	}
	
	public String toString(){
		String retval="multiplicity: ";
		for(PopBin b:mmap.keySet())
			retval+="<"+b.getId() + ":" + mmap.get(b)+"> ";
		return retval;
	}
	public boolean isEmpty() {
		if(mmap.isEmpty()||getSum()==0)
			return true;
		else 
			return false;
	}
}
