package org.clas.analysis;

//  @author m.c.kunkel
import java.util.Arrays;
import java.util.HashMap;
import java.util.Map;

import org.jlab.groot.data.H1F;
import org.jlab.groot.math.F1D;
import org.jlab.groot.math.RandomFunc;
import org.jlab.groot.ui.TCanvas;

public class Coordinate {
	private Integer[] size;

	public Coordinate(Integer... size) {
		this.size = size;
	}

	@Override
	public int hashCode() {
		final int prime = 31;
		int result = 1;
		result = prime * result + Arrays.hashCode(size);
		return result;
	}

	@Override
	public boolean equals(Object obj) {
		if (this == obj)
			return true;
		if (obj == null)
			return false;
		if (getClass() != obj.getClass())
			return false;
		Coordinate other = (Coordinate) obj;
		if (!Arrays.equals(size, other.size))
			return false;
		return true;
	}

	public static void main(String[] args) {
		Map<Coordinate, H1F> trial = new HashMap<Coordinate, H1F>();

		trial.put(new Coordinate(1, 1), new H1F("h" + 1, "", 50, 0, 2.0 * Math.PI));
		trial.put(new Coordinate(2, 1), new H1F("h" + 2, "", 50, 0, 2.0 * Math.PI));
		trial.put(new Coordinate(1, 1, 1, 1, 1, 1), new H1F("h1111" + 1, "", 50, 0, 2.0 * Math.PI));
		trial.put(new Coordinate(1, 1, 1), new H1F("h11" + 1, "", 50, 0, 2.0 * Math.PI));

		F1D func = new F1D("func", "2.0+[a]*cos(x)+[b]*cos(2*x)", 0.0, 2.0 * 3.1415);
		F1D func2 = new F1D("func2", "2.0+[a]*sin(x)+[b]*sin(2*x)", 0.0, 2.0 * 3.1415);
		F1D func3 = new F1D("func3", "2.0+[a]*sin(x)", 0.0, 2.0 * 3.1415);
		F1D func4 = new F1D("func4", "2.0+[a]*cos(x)", 0.0, 2.0 * 3.1415);

		func.setParameter(0, 0.5);
		func.setParameter(1, 1.0);
		func2.setParameter(0, 0.5);
		func2.setParameter(1, 1.0);
		func3.setParameter(0, 0.5);
		func4.setParameter(0, 0.5);
		RandomFunc randfunc = new RandomFunc(func);
		RandomFunc randfunc2 = new RandomFunc(func2);
		RandomFunc randfunc3 = new RandomFunc(func3);
		RandomFunc randfunc4 = new RandomFunc(func4);

		for (int i = 0; i < 1800; i++) {
			trial.get(new Coordinate(1, 1)).fill(randfunc.random());
			trial.get(new Coordinate(2, 1)).fill(randfunc2.random());
			trial.get(new Coordinate(1, 1, 1, 1, 1, 1)).fill(randfunc3.random());
			trial.get(new Coordinate(1, 1, 1)).fill(randfunc4.random());
		}

		trial.get(new Coordinate(1, 1)).setTitleX("First");
		trial.get(new Coordinate(2, 1)).setTitleX("Second");
		trial.get(new Coordinate(1, 1, 1, 1, 1, 1)).setTitleX("Third");
                trial.get(new Coordinate(1, 1, 1)).setTitleX("Fourth");

		TCanvas c1 = new TCanvas("c1", 800, 800);
		c1.divide(2, 2);
		c1.cd(0);
		c1.draw(trial.get(new Coordinate(1, 1)));
		c1.cd(1);
		c1.draw(trial.get(new Coordinate(2, 1)));
		c1.cd(2);
		c1.draw(trial.get(new Coordinate(1, 1, 1, 1, 1, 1)));
		c1.cd(3);
		c1.draw(trial.get(new Coordinate(1, 1, 1)));

	}

}