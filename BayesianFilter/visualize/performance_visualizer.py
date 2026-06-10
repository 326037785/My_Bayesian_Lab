"""
性能可视化模块
"""
import numpy as np
import matplotlib.pyplot as plt
from typing import Optional, List, Dict, Tuple


class PerformanceVisualizer:
    """性能可视化器
    
    用于可视化跟踪性能指标，包括：
    - RMSE时序图
    - OSPA/GOSPA对比图
    - 误差分布图
    - 性能指标统计
    """
    
    def __init__(self, 
                 figsize: Tuple[int, int] = (12, 8),
                 style: str = 'seaborn-v0_8-darkgrid'):
        """
        初始化性能可视化器
        
        Args:
            figsize: 图形大小
            style: matplotlib样式
        """
        self.figsize = figsize
        self.style = style
        
        # 颜色配置
        self.colors = plt.cm.tab10(np.linspace(0, 1, 10))
        
        # 设置样式
        try:
            plt.style.use(style)
        except:
            plt.style.use('default')
    
    def plot_rmse_comparison(self,
                             rmse_dict: Dict[str, List[float]],
                             time_steps: Optional[np.ndarray] = None,
                             title: str = "RMSE Comparison",
                             xlabel: str = "Time Step",
                             ylabel: str = "RMSE",
                             save_path: Optional[str] = None) -> plt.Figure:
        """绘制RMSE对比图
        
        Args:
            rmse_dict: RMSE字典，键为算法名称，值为RMSE列表
            time_steps: 时间步数组
            title: 图形标题
            xlabel: x轴标签
            ylabel: y轴标签
            save_path: 保存路径
            
        Returns:
            matplotlib图形对象
        """
        fig, ax = plt.subplots(figsize=self.figsize)
        
        for i, (name, rmse_values) in enumerate(rmse_dict.items()):
            color_idx = i % len(self.colors)
            
            if time_steps is not None:
                x = time_steps[:len(rmse_values)]
            else:
                x = np.arange(len(rmse_values))
            
            ax.plot(x, rmse_values, color=self.colors[color_idx], 
                   linewidth=2, label=name, alpha=0.8)
        
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        ax.set_title(title)
        ax.legend(loc='best')
        ax.grid(True, alpha=0.3)
        
        plt.tight_layout()
        
        if save_path:
            plt.savefig(save_path, dpi=300, bbox_inches='tight')
        
        return fig
    
    def plot_ospa_comparison(self,
                              ospa_dict: Dict[str, Tuple[List[float], List[float]]],
                              time_steps: Optional[np.ndarray] = None,
                              title: str = "OSPA Comparison",
                              save_path: Optional[str] = None) -> plt.Figure:
        """绘制OSPA对比图
        
        Args:
            ospa_dict: OSPA字典，键为算法名称，值为(总OSPA, 基数误差)元组
            time_steps: 时间步数组
            title: 图形标题
            save_path: 保存路径
            
        Returns:
            matplotlib图形对象
        """
        fig, axes = plt.subplots(2, 1, figsize=self.figsize, sharex=True)
        
        for i, (name, (ospa_values, card_errors)) in enumerate(ospa_dict.items()):
            color_idx = i % len(self.colors)
            
            if time_steps is not None:
                x = time_steps[:len(ospa_values)]
            else:
                x = np.arange(len(ospa_values))
            
            axes[0].plot(x, ospa_values, color=self.colors[color_idx], 
                        linewidth=2, label=name, alpha=0.8)
            axes[1].plot(x, card_errors, color=self.colors[color_idx], 
                        linewidth=2, label=name, alpha=0.8)
        
        axes[0].set_ylabel('OSPA')
        axes[0].set_title(title)
        axes[0].legend(loc='best')
        axes[0].grid(True, alpha=0.3)
        
        axes[1].set_xlabel('Time Step')
        axes[1].set_ylabel('Cardinality Error')
        axes[1].legend(loc='best')
        axes[1].grid(True, alpha=0.3)
        
        plt.tight_layout()
        
        if save_path:
            plt.savefig(save_path, dpi=300, bbox_inches='tight')
        
        return fig
    
    def plot_gospa_comparison(self,
                               gospa_dict: Dict[str, Dict[str, List[float]]],
                               time_steps: Optional[np.ndarray] = None,
                               title: str = "GOSPA Comparison",
                               save_path: Optional[str] = None) -> plt.Figure:
        """绘制GOSPA对比图
        
        Args:
            gospa_dict: GOSPA字典，键为算法名称，值为包含'gospa', 'localization', 'missed', 'false_alarm'的字典
            time_steps: 时间步数组
            title: 图形标题
            save_path: 保存路径
            
        Returns:
            matplotlib图形对象
        """
        fig, axes = plt.subplots(2, 2, figsize=self.figsize)
        axes = axes.flatten()
        
        metrics = ['gospa', 'localization', 'missed', 'false_alarm']
        metric_titles = ['Total GOSPA', 'Localization Error', 'Missed Detection', 'False Alarm']
        
        for i, (name, metric_dict) in enumerate(gospa_dict.items()):
            color_idx = i % len(self.colors)
            
            for j, metric in enumerate(metrics):
                if metric in metric_dict:
                    values = metric_dict[metric]
                    
                    if time_steps is not None:
                        x = time_steps[:len(values)]
                    else:
                        x = np.arange(len(values))
                    
                    axes[j].plot(x, values, color=self.colors[color_idx], 
                               linewidth=2, label=name, alpha=0.8)
        
        for j, (metric, metric_title) in enumerate(zip(metrics, metric_titles)):
            axes[j].set_xlabel('Time Step')
            axes[j].set_ylabel(metric_title)
            axes[j].set_title(metric_title)
            axes[j].legend(loc='best')
            axes[j].grid(True, alpha=0.3)
        
        fig.suptitle(title, fontsize=14)
        plt.tight_layout()
        
        if save_path:
            plt.savefig(save_path, dpi=300, bbox_inches='tight')
        
        return fig
    
    def plot_error_distribution(self,
                                 errors: Dict[str, np.ndarray],
                                 title: str = "Error Distribution",
                                 xlabel: str = "Error",
                                 ylabel: str = "Frequency",
                                 bins: int = 30,
                                 save_path: Optional[str] = None) -> plt.Figure:
        """绘制误差分布图
        
        Args:
            errors: 误差字典，键为算法名称，值为误差数组
            title: 图形标题
            xlabel: x轴标签
            ylabel: y轴标签
            bins: 直方图bin数量
            save_path: 保存路径
            
        Returns:
            matplotlib图形对象
        """
        fig, ax = plt.subplots(figsize=self.figsize)
        
        for i, (name, error_values) in enumerate(errors.items()):
            color_idx = i % len(self.colors)
            ax.hist(error_values, bins=bins, alpha=0.5, 
                   color=self.colors[color_idx], label=name, density=True)
        
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        ax.set_title(title)
        ax.legend(loc='best')
        ax.grid(True, alpha=0.3)
        
        plt.tight_layout()
        
        if save_path:
            plt.savefig(save_path, dpi=300, bbox_inches='tight')
        
        return fig
    
    def plot_performance_summary(self,
                                  metrics_dict: Dict[str, Dict[str, float]],
                                  title: str = "Performance Summary",
                                  save_path: Optional[str] = None) -> plt.Figure:
        """绘制性能摘要图
        
        Args:
            metrics_dict: 指标字典，键为算法名称，值为指标字典
            title: 图形标题
            save_path: 保存路径
            
        Returns:
            matplotlib图形对象
        """
        fig, axes = plt.subplots(2, 2, figsize=self.figsize)
        
        # 提取指标
        algorithms = list(metrics_dict.keys())
        metric_names = ['rmse_mean', 'ospa_mean', 'gospa_mean', 'cardinality_error_mean']
        metric_labels = ['RMSE', 'OSPA', 'GOSPA', 'Cardinality Error']
        
        for j, (metric_name, metric_label) in enumerate(zip(metric_names, metric_labels)):
            values = []
            for alg in algorithms:
                if metric_name in metrics_dict[alg]:
                    values.append(metrics_dict[alg][metric_name])
                else:
                    values.append(0)
            
            x = np.arange(len(algorithms))
            bars = axes[j // 2, j % 2].bar(x, values, alpha=0.7)
            
            # 为每个bar设置不同颜色
            for k, bar in enumerate(bars):
                bar.set_color(self.colors[k % len(self.colors)])
            
            axes[j // 2, j % 2].set_xlabel('Algorithm')
            axes[j // 2, j % 2].set_ylabel(metric_label)
            axes[j // 2, j % 2].set_title(metric_label)
            axes[j // 2, j % 2].set_xticks(x)
            axes[j // 2, j % 2].set_xticklabels(algorithms, rotation=45, ha='right')
            axes[j // 2, j % 2].grid(True, alpha=0.3)
        
        fig.suptitle(title, fontsize=14)
        plt.tight_layout()
        
        if save_path:
            plt.savefig(save_path, dpi=300, bbox_inches='tight')
        
        return fig
    
    def plot_cardinality_over_time(self,
                                    true_cardinality: List[int],
                                    estimated_cardinality: Dict[str, List[int]],
                                    time_steps: Optional[np.ndarray] = None,
                                    title: str = "Cardinality Over Time",
                                    save_path: Optional[str] = None) -> plt.Figure:
        """绘制目标数量随时间变化图
        
        Args:
            true_cardinality: 真实目标数量列表
            estimated_cardinality: 估计目标数量字典
            time_steps: 时间步数组
            title: 图形标题
            save_path: 保存路径
            
        Returns:
            matplotlib图形对象
        """
        fig, ax = plt.subplots(figsize=self.figsize)
        
        if time_steps is not None:
            x = time_steps[:len(true_cardinality)]
        else:
            x = np.arange(len(true_cardinality))
        
        ax.plot(x, true_cardinality, 'k-', linewidth=2, label='True', alpha=0.8)
        
        for i, (name, est_card) in enumerate(estimated_cardinality.items()):
            color_idx = i % len(self.colors)
            ax.plot(x[:len(est_card)], est_card, '--', 
                   color=self.colors[color_idx], linewidth=1.5, 
                   label=name, alpha=0.7)
        
        ax.set_xlabel('Time Step')
        ax.set_ylabel('Number of Targets')
        ax.set_title(title)
        ax.legend(loc='best')
        ax.grid(True, alpha=0.3)
        
        plt.tight_layout()
        
        if save_path:
            plt.savefig(save_path, dpi=300, bbox_inches='tight')
        
        return fig
    
    def plot_computation_time(self,
                               time_dict: Dict[str, float],
                               title: str = "Computation Time",
                               save_path: Optional[str] = None) -> plt.Figure:
        """绘制计算时间对比图
        
        Args:
            time_dict: 时间字典，键为算法名称，值为计算时间（秒）
            title: 图形标题
            save_path: 保存路径
            
        Returns:
            matplotlib图形对象
        """
        fig, ax = plt.subplots(figsize=self.figsize)
        
        algorithms = list(time_dict.keys())
        times = list(time_dict.values())
        
        x = np.arange(len(algorithms))
        bars = ax.bar(x, times, alpha=0.7)
        
        for i, bar in enumerate(bars):
            bar.set_color(self.colors[i % len(self.colors)])
        
        ax.set_xlabel('Algorithm')
        ax.set_ylabel('Time (seconds)')
        ax.set_title(title)
        ax.set_xticks(x)
        ax.set_xticklabels(algorithms, rotation=45, ha='right')
        ax.grid(True, alpha=0.3)
        
        plt.tight_layout()
        
        if save_path:
            plt.savefig(save_path, dpi=300, bbox_inches='tight')
        
        return fig
