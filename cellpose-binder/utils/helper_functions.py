def save_results(masks, image_path, output_dir='./results'):
    """Save segmentation results"""
    import os
    import numpy as np
    from PIL import Image
    
    os.makedirs(output_dir, exist_ok=True)
    
    # Save masks
    mask_path = os.path.join(output_dir, 'masks.png')
    Image.fromarray(masks.astype(np.uint16)).save(mask_path)
    
    print(f"💾 Results saved to {output_dir}")
    return mask_path

def analyze_results(masks):
    """Simple analysis of segmentation results"""
    import numpy as np
    
    unique_masks = np.unique(masks)
    n_cells = len(unique_masks) - 1  # exclude background
    
    cell_sizes = []
    for cell_id in unique_masks[1:]:  # skip background
        cell_size = np.sum(masks == cell_id)
        cell_sizes.append(cell_size)
    
    print(f"📊 Analysis Results:")
    print(f"   Number of cells: {n_cells}")
    print(f"   Average cell size: {np.mean(cell_sizes):.1f} pixels")
    print(f"   Size range: {np.min(cell_sizes)} - {np.max(cell_sizes)} pixels")
    
    return {
        'n_cells': n_cells,
        'cell_sizes': cell_sizes,
        'avg_size': np.mean(cell_sizes)
    }