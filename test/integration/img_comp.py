import numpy as np
import matplotlib.pyplot as plt

def img_to_tiles(img):
    M = (img.shape[0] // 4) + 1
    N = (img.shape[1] // 4) + 1
    tiles = [img[x:x+M,y:y+N] for x in range(0,img.shape[0],M) for y in range(0,img.shape[1],N)]
    return tiles

def img_to_hist(img):
    bins = np.arange(-.1, 1.1, .1)
    hist = np.array([])
    for channel in range(3):
        hist = np.append(
            hist,
            np.histogram(img[channel], bins=bins)[0] / img[channel].size
        )
    return hist

def tiles_to_hists(tiles):
    hists = np.array([])
    for tile in tiles:
        hists = np.append(hists, img_to_hist(tile))
    return hists

def hists_approx_equal(hist1, hist2, eps=.2):
    for freq1, freq2 in zip(hist1, hist2):
        if np.abs(freq1 - freq2) > eps:
            return False
    return True

def imgs_approx_equal(img1, img2):
    if isinstance(img1, str):
        img1 = plt.imread(img1)
    if isinstance(img2, str):
        img2 = plt.imread(img2)
    if img1.shape != img2.shape:
        return False
    return hists_approx_equal(
        tiles_to_hists(img_to_tiles(img1)),
        tiles_to_hists(img_to_tiles(img2)),
    )
