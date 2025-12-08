import numpy as np

class AuxiliaryStatistics:
    """
    Computes comprehensive auxiliary statistics for the test.
    """
    
    @staticmethod
    def _compute_all(data: np.ndarray) -> np.ndarray:
        """Internal method to compute all statistics once."""
        # Base percentiles
        percentiles = np.percentile(data, [25, 50, 75, 90, 95])
        q_ratios = [
            percentiles[2] / percentiles[0], 
            percentiles[3] / percentiles[1],
            percentiles[3] / percentiles[0], 
            percentiles[4] / percentiles[2]
        ]

        # Additional percentile sets
        percentiles2 = np.percentile(data, [5,10,15,20,30,35,40,45,55,60,65,70,80,85])
        percentiles3 = np.percentile(data, [4,8,12,16,24,28,32,36,44,48,52,56,64,68,72])
        percentiles4 = np.percentile(data, [7,14,21,42,49,63,77,84,91,98])
        percentiles5 = np.percentile(data, [9,19,29,39,49,59,69,79,89,99])

        # Quantile ratios
        q_ratio2 = [
            percentiles2[5] / percentiles2[1], percentiles2[9] / percentiles2[3],
            percentiles2[9] / percentiles2[1], percentiles2[11] / percentiles2[9],
            percentiles2[11] / percentiles2[3], percentiles2[11] / percentiles2[1],
            percentiles2[13] / percentiles2[11], percentiles2[13] / percentiles2[9],
            percentiles2[13] / percentiles2[3], percentiles2[13] / percentiles2[1],
            percentiles2[13] / percentiles2[5], percentiles2[13] / percentiles2[7],
            percentiles2[13] / percentiles2[4], percentiles2[13] / percentiles2[0],
            percentiles2[13] / percentiles2[2]
        ]

        q_ratio3 = [
            percentiles3[4] / percentiles3[0], percentiles3[6] / percentiles3[2],
            percentiles3[8] / percentiles3[4], percentiles3[10] / percentiles3[6],
            percentiles3[12] / percentiles3[8], percentiles3[14] / percentiles3[10],
            percentiles3[14] / percentiles3[0], percentiles3[13] / percentiles3[1],
            percentiles3[11] / percentiles3[3], percentiles3[7] / percentiles3[5]
        ]

        q_ratio4 = [
            percentiles4[3] / percentiles4[0], percentiles4[4] / percentiles4[1],
            percentiles4[5] / percentiles4[2], percentiles4[6] / percentiles4[3],
            percentiles4[7] / percentiles4[4], percentiles4[8] / percentiles4[5],
            percentiles4[9] / percentiles4[6], percentiles4[9] / percentiles4[0],
            percentiles4[8] / percentiles4[2], percentiles4[7] / percentiles4[1]
        ]

        q_ratio5 = [
            percentiles5[2] / percentiles5[0], percentiles5[3] / percentiles5[1],
            percentiles5[4] / percentiles5[2], percentiles5[5] / percentiles5[3],
            percentiles5[5] / percentiles5[0], percentiles5[4] / percentiles5[1],
            percentiles5[3] / percentiles5[0]
        ]

        # Mixed and combined ratios
        q_ratio_mixed = [
            percentiles5[9] / percentiles4[0], percentiles5[8] / percentiles4[1],
            percentiles4[9] / percentiles5[0], percentiles4[8] / percentiles5[1],
            percentiles5[5] / percentiles4[5], percentiles4[4] / percentiles5[4],
            percentiles5[2] / percentiles4[2], percentiles4[3] / percentiles5[3],
            percentiles5[6] / percentiles4[6], percentiles4[7] / percentiles5[7]
        ]

        q_ratio_combined = [
            percentiles3[14] / percentiles2[0], percentiles2[13] / percentiles3[0],
            percentiles2[12] / percentiles3[1], percentiles3[13] / percentiles2[2],
            percentiles3[12] / percentiles2[3], percentiles2[11] / percentiles3[3],
            percentiles3[11] / percentiles2[4], percentiles2[9] / percentiles3[6],
            percentiles3[8] / percentiles2[7], percentiles3[4] / percentiles2[5]
        ]

        # Inequality measures
        # Gini coefficient
        # Mean absolute difference
        mad = np.abs(np.subtract.outer(data, data)).mean()
        gini = mad / (2 * np.mean(data))
        
        # Theil index
        mean_data = np.mean(data)
        theil = np.mean((data / mean_data) * np.log(data / mean_data))

        # Fine-grained percentiles
        fine_percentiles = np.percentile(data, np.arange(1, 100))
        fine_ratios = fine_percentiles[1:] / fine_percentiles[:-1]

        # Combine all statistics
        stats = np.hstack([
            percentiles, q_ratios, gini, theil,
            percentiles2, q_ratio2,
            percentiles3, q_ratio3, q_ratio_combined,
            percentiles4, q_ratio4,
            percentiles5, q_ratio5,
            q_ratio_mixed,
            fine_percentiles, fine_ratios
        ])
        return stats

    @staticmethod
    def compute(data: np.ndarray, order: int) -> float:
        """
        Compute the auxiliary statistic for a given order.
        
        Args:
            data (np.ndarray): Input data.
            order (int): The index of the statistic to return (1-based).
            
        Returns:
            float: The computed statistic.
        """
        stats = AuxiliaryStatistics._compute_all(data)

        if order < 1 or order > len(stats):
            raise ValueError(f"Order must be between 1 and {len(stats)} (you passed {order})")
        
        return stats[order - 1]
    
    @staticmethod
    def compute_vector(data: np.ndarray, max_order: int) -> np.ndarray:
        """
        Compute all statistics up to max_order efficiently.
        """
        stats = AuxiliaryStatistics._compute_all(data)
        
        if max_order > len(stats):
             raise ValueError(f"Max order {max_order} exceeds available statistics {len(stats)}")
             
        return stats[:max_order]
