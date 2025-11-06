package gedi.iTiSS.utils.multithreading;

import java.util.HashSet;
import java.util.Set;

public abstract class NotifyOnFinishedRunnable implements Runnable {
    private Set<RunnableFinishedListener> listeners = new HashSet<>();

    public final void addListener(RunnableFinishedListener listener) {
        this.listeners.add(listener);
    }

    public final void removeListener(RunnableFinishedListener listener) {
        this.listeners.remove(listener);
    }

    private void notifyListeners() {
        for (RunnableFinishedListener listener : listeners) {
            listener.notifyRunnableFinished(this);
        }
    }

    @Override
    public final void run() {
//        try {
//            doRun();
//        } finally {
//            notifyListeners();
//        }
        doRun();
        notifyListeners();
    }

    public abstract void doRun();
}
