import numpy as np

def half_max_positions(ap_bins, avg_mrna):
    ap_bins = np.array(ap_bins)
    avg_mrna = np.array(avg_mrna)
    i_max = np.nanargmax(avg_mrna)
    max_x = ap_bins[i_max]
    max_y = avg_mrna[i_max]
    half = max_y / 2

    # Left side
    left_idx = np.where(ap_bins < max_x)[0]
    left_half = np.nan
    if left_idx.size > 0:
        # find last point above half on the left
        mask = avg_mrna[left_idx] >= half
        if np.any(mask):
            j = left_idx[np.where(mask)[0][-1]]
            if j+1 <= i_max-1:  # make sure there is a point below half to interpolate with
                x1, y1 = ap_bins[j], avg_mrna[j]
                x2, y2 = ap_bins[j+1], avg_mrna[j+1]
                # linear interpolation to half
                left_half = x1 + (half - y1) * (x2 - x1) / (y2 - y1)

    # Right side
    right_idx = np.where(ap_bins > max_x)[0]
    right_half = np.nan
    if right_idx.size > 0:
        mask = avg_mrna[right_idx] >= half
        if np.any(mask):
            j = right_idx[np.where(mask)[0][0]] - 1 + 1  # first above-half to the right (careful indexing)
            if j-1 >= i_max+1:
                x1, y1 = ap_bins[j-1], avg_mrna[j-1]
                x2, y2 = ap_bins[j], avg_mrna[j]
                right_half = x1 + (half - y1) * (x2 - x1) / (y2 - y1)

    return max_x, max_y, left_half, right_half, half
