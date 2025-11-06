package gedi.iTiSS.utils.loader;

import gedi.core.region.GenomicRegionStoragePreload;
import gedi.core.region.intervalTree.MemoryIntervalTreeStorage;
import gedi.core.workspace.loader.WorkspaceItemLoader;
import gedi.util.dynamic.DynamicObject;

import java.io.IOException;
import java.nio.file.Path;
import java.util.HashSet;

public class TsrFileLoader implements WorkspaceItemLoader<MemoryIntervalTreeStorage<TsrData>, GenomicRegionStoragePreload> {

    public static final String[] extensions = new String[]{"tsr"};

    @Override
    public String[] getExtensions() {
        return new String[0];
    }

    @Override
    public MemoryIntervalTreeStorage<TsrData> load(Path path) throws IOException {
        return new TsrDataFileReader(path.toString()).readIntoMemoryTakeFirst();
    }

    @Override
    public GenomicRegionStoragePreload preload(Path path) throws IOException {
        return new GenomicRegionStoragePreload(TsrData.class, new TsrData(0, 10, new HashSet<>()), DynamicObject.getEmpty());
    }

    @Override
    public Class<MemoryIntervalTreeStorage<TsrData>> getItemClass() {
        return (Class)MemoryIntervalTreeStorage.class;
    }

    @Override
    public boolean hasOptions() {
        return false;
    }

    @Override
    public void updateOptions(Path path) {

    }
}
