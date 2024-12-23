package gedi.util.math.stat.counting;

import java.util.Arrays;

public class ItemCount<T> {

	private T item;
	private int[] count;
	public ItemCount(T item, int[] count) {
		this.item = item;
		this.count = count;
	}
	public T getItem() {
		return item;
	}
	public int[] getCount() {
		return count;
	}
	
	public int getCount(int d) {
		return count[d];
	}
	
	public int getDimensions() {
		return count.length;
	}
	@Override
	public String toString() {
		return "ItemCount [item=" + item + ", count=" + Arrays.toString(count)
				+ "]";
	}
	@Override
	public int hashCode() {
		final int prime = 31;
		int result = 1;
		result = prime * result + Arrays.hashCode(count);
		result = prime * result + ((item == null) ? 0 : item.hashCode());
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
		ItemCount other = (ItemCount) obj;
		if (!Arrays.equals(count, other.count))
			return false;
		if (item == null) {
			if (other.item != null)
				return false;
		} else if (!item.equals(other.item))
			return false;
		return true;
	}

	
	
	
}
