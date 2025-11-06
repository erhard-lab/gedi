package gedi.iTiSS.utils.scheduling;

import java.util.Arrays;

public class AccessionList {
    private AccessionStatus[] accession;
    private int nextNeeded;

    public AccessionList(int accessionSize) {
        accession = new AccessionStatus[accessionSize];
        Arrays.fill(accession, AccessionStatus.NOT_ACCESSED);
        nextNeeded = 0;
    }

//    public int getNextNeeded() {
////        if (accession[nextNeeded] == AccessionStatus.NOT_ACCESSED) {
////            return nextNeeded;
////        }
//        int tmp = nextNeeded+1;
//        while (tmp != nextNeeded) {
//            if (tmp == accession.length) {
//                tmp = 0;
//            }
//            if (accession[tmp] == AccessionStatus.NOT_ACCESSED) {
//                nextNeeded = tmp;
//                return nextNeeded;
//            }
//            tmp++;
//        }
//        return -1;
//    }

    public int getNextNeeded() {
        for (int i = nextNeeded+1; i < accession.length; i++) {
            if (accession[i] == AccessionStatus.NOT_ACCESSED) {
                nextNeeded = i;
                return nextNeeded;
            }
        }
        return -1;
    }

    public int getNextNeededReset() {
        for (int i = 0; i < accession.length; i++) {
            if (accession[i] == AccessionStatus.NOT_ACCESSED) {
                nextNeeded = i;
                return nextNeeded;
            }
        }
        return -1;
    }

    public void stopAccess(int index) {
        accession[index] = AccessionStatus.FINISHED;
    }

    public void startAccess(int index) {
        accession[index] = AccessionStatus.CURRENTLY_ACCESSED;
    }

    public boolean isFree(int index) {
        return accession[index] != AccessionStatus.CURRENTLY_ACCESSED;
    }

    public boolean allFinished() {
        for (AccessionStatus anAccession : accession) {
            if (anAccession != AccessionStatus.FINISHED) {
                return false;
            }
        }
        return true;
    }
}
