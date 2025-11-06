package gedi.iTiSS.data;

import java.util.ArrayList;
import java.util.List;

/**
 * This class is used as storage class for accurate CIT-access regarding the ordering.
 * It contains two access lists with the indices of the CIT file in the first {@code citAccess} list and the
 * lanes that should be accessed in that CIT-file in the second list {@code laneAccess}.
 * Exp.:
 * CIT-files: 0, 1, 2 (each 2 lanes)
 * The files should be accessed in the following order: 2.0, 1.1, 2.1, 0.0, 0.1
 * The {@code citAccess} list of this class then would contain: [2, 1, 2, 0]
 * The {@code laneAccess} list of this class then would contain: [[0], [1], [1], [0,1]]
 */
public class CitAccessInfo {
    private List<Integer> citAccess;
    private List<List<Integer>> laneAccess;

    public CitAccessInfo() {
        this.citAccess = new ArrayList<>();
        this.laneAccess = new ArrayList<>();
    }

    public void add(int cit, int lane) {
        if (citAccess.size() == 0 || citAccess.get(citAccess.size()-1) != cit) {
            this.citAccess.add(cit);
            List<Integer> l = new ArrayList<>();
            l.add(lane);
            laneAccess.add(l);
        }
        else {
            laneAccess.get(citAccess.size()-1).add(lane);
        }
    }

    public int getCitAccess(int index) {
        return citAccess.get(index);
    }

    public List<Integer> getLaneAccess(int cit) {
        return laneAccess.get(cit);
    }

    public int getCitAccessNum() {
        return citAccess.size();
    }
}
