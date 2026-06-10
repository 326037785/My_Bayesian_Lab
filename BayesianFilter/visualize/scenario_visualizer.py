"""
场景可视化模块
"""
import numpy as np
import matplotlib.pyplot as plt
from typing import Optional, List, Dict, Tuple
from matplotlib.patches import Ellipse
import matplotlib.animation as animation


class ScenarioVisualizer:
    """场景可视化器
    
    用于可视化跟踪场景，包括：
    - 真实轨迹
    - 观测点
    - 滤波结果
    - 不确定性椭圆
    """
    
    def __init__(self, 
                 figsize: Tuple[int, int] = (10, 8),
                 style: str = 'seaborn-v0_8-darkgrid'):
        """
        初始化场景可视化器
        
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
    
    def plot_scenario(self,
                      true_trajectories: Optional[Dict[int, np.ndarray]] = None,
                      measurements: Optional[List[np.ndarray]] = None,
                      estimated_trajectories: Optional[Dict[int, np.ndarray]] = None,
                      title: str = "Tracking Scenario",
                      show_measurements: bool = True,
                      show_estimates: bool = True,
                      show_true: bool = True,
                      uncertainty_ellipses: Optional[Dict[int, List[np.ndarray]]] = None,
                      save_path: Optional[str] = None) -> plt.Figure:
        """绘制跟踪场景
        
        Args:
            true_trajectories: 真实轨迹字典，键为目标ID，值为轨迹数组
            measurements: 观测列表，每个元素为一个时刻的观测数组
            estimated_trajectories: 估计轨迹字典
            title: 图形标题
            show_measurements: 是否显示观测
            show_estimates: 是否显示估计
            show_true: 是否显示真实轨迹
            uncertainty_ellipses: 不确定性椭圆字典
            save_path: 保存路径
            
        Returns:
            matplotlib图形对象
        """
        fig, ax = plt.subplots(figsize=self.figsize)
        
        # 绘制真实轨迹
        if show_true and true_trajectories is not None:
            for target_id, trajectory in true_trajectories.items():
                color_idx = target_id % len(self.colors)
                ax.plot(trajectory[:, 0], trajectory[:, 1], 
                       color=self.colors[color_idx], linewidth=2, 
                       label=f'True {target_id}', alpha=0.8)
                
                # 标记起点和终点
                ax.plot(trajectory[0, 0], trajectory[0, 1], 
                       'o', color=self.colors[color_idx], markersize=8)
                ax.plot(trajectory[-1, 0], trajectory[-1, 1], 
                       's', color=self.colors[color_idx], markersize=8)
        
        # 绘制观测
        if show_measurements and measurements is not None:
            all_meas = []
            for meas in measurements:
                if meas is not None and len(meas) > 0:
                    all_meas.append(meas)
            
            if all_meas:
                all_meas = np.vstack(all_meas)
                ax.scatter(all_meas[:, 0], all_meas[:, 1], 
                          c='gray', marker='x', s=20, alpha=0.3, 
                          label='Measurements')
        
        # 绘制估计轨迹
        if show_estimates and estimated_trajectories is not None:
            for target_id, trajectory in estimated_trajectories.items():
                color_idx = target_id % len(self.colors)
                ax.plot(trajectory[:, 0], trajectory[:, 1], 
                       '--', color=self.colors[color_idx], linewidth=1.5,
                       label=f'Estimated {target_id}', alpha=0.7)
        
        # 绘制不确定性椭圆
        if uncertainty_ellipses is not None:
            for target_id, ellipses in uncertainty_ellipses.items():
                color_idx = target_id % len(self.colors)
                for ellipse_params in ellipses:
                    # ellipse_params: [x, y, width, height, angle]
                    ellipse = Ellipse(
                        (ellipse_params[0], ellipse_params[1]),
                        ellipse_params[2], ellipse_params[3],
                        angle=np.degrees(ellipse_params[4]),
                        fill=False, 
                        edgecolor=self.colors[color_idx],
                        linestyle='--', alpha=0.3
                    )
                    ax.add_patch(ellipse)
        
        # 设置图形属性
        ax.set_xlabel('X Position')
        ax.set_ylabel('Y Position')
        ax.set_title(title)
        
        # 当目标数量过多时，不显示图例（避免覆盖图形）
        n_true = len(true_trajectories) if true_trajectories else 0
        n_est = len(estimated_trajectories) if estimated_trajectories else 0
        if n_true + n_est <= 20:
            ax.legend(loc='best', fontsize='small')
        
        ax.grid(True, alpha=0.3)
        ax.set_aspect('equal')
        
        plt.tight_layout()
        
        if save_path:
            plt.savefig(save_path, dpi=300, bbox_inches='tight')
        
        return fig
    
    def plot_measurements_only(self,
                                measurements: List[np.ndarray],
                                title: str = "Measurements",
                                save_path: Optional[str] = None) -> plt.Figure:
        """仅绘制观测
        
        Args:
            measurements: 观测列表
            title: 图形标题
            save_path: 保存路径
            
        Returns:
            matplotlib图形对象
        """
        fig, ax = plt.subplots(figsize=self.figsize)
        
        all_meas = []
        for t, meas in enumerate(measurements):
            if meas is not None and len(meas) > 0:
                for m in meas:
                    all_meas.append(m)
        
        if all_meas:
            all_meas = np.array(all_meas)
            ax.scatter(all_meas[:, 0], all_meas[:, 1], 
                      c='blue', marker='x', s=30, alpha=0.5)
        
        ax.set_xlabel('X Position')
        ax.set_ylabel('Y Position')
        ax.set_title(title)
        ax.grid(True, alpha=0.3)
        ax.set_aspect('equal')
        
        plt.tight_layout()
        
        if save_path:
            plt.savefig(save_path, dpi=300, bbox_inches='tight')
        
        return fig
    
    def animate_scenario(self,
                          true_trajectories: Dict[int, np.ndarray],
                          measurements: List[np.ndarray],
                          estimated_trajectories: Optional[Dict[int, np.ndarray]] = None,
                          interval: int = 200,
                          save_path: Optional[str] = None) -> animation.FuncAnimation:
        """创建跟踪场景动画
        
        Args:
            true_trajectories: 真实轨迹字典
            measurements: 观测列表
            estimated_trajectories: 估计轨迹字典
            interval: 帧间隔（毫秒）
            save_path: 保存路径
            
        Returns:
            matplotlib动画对象
        """
        fig, ax = plt.subplots(figsize=self.figsize)
        
        # 获取时间步数
        n_steps = len(measurements)
        
        # 计算坐标范围
        all_positions = []
        for trajectory in true_trajectories.values():
            all_positions.append(trajectory)
        if estimated_trajectories:
            for trajectory in estimated_trajectories.values():
                all_positions.append(trajectory)
        
        if all_positions:
            all_positions = np.vstack(all_positions)
            x_min, x_max = all_positions[:, 0].min() - 50, all_positions[:, 0].max() + 50
            y_min, y_max = all_positions[:, 1].min() - 50, all_positions[:, 1].max() + 50
        else:
            x_min, x_max = -500, 500
            y_min, y_max = -500, 500
        
        # 初始化图形元素
        true_lines = {}
        est_lines = {}
        meas_scatter = ax.scatter([], [], c='gray', marker='x', s=20, alpha=0.5)
        
        for target_id in true_trajectories.keys():
            color_idx = target_id % len(self.colors)
            line, = ax.plot([], [], color=self.colors[color_idx], linewidth=2, 
                          label=f'True {target_id}')
            true_lines[target_id] = line
            
            if estimated_trajectories and target_id in estimated_trajectories:
                line, = ax.plot([], [], '--', color=self.colors[color_idx], 
                              linewidth=1.5, label=f'Est {target_id}')
                est_lines[target_id] = line
        
        ax.set_xlim(x_min, x_max)
        ax.set_ylim(y_min, y_max)
        ax.set_xlabel('X Position')
        ax.set_ylabel('Y Position')
        ax.set_title('Tracking Animation')
        ax.legend(loc='best')
        ax.grid(True, alpha=0.3)
        ax.set_aspect('equal')
        
        def init():
            for line in true_lines.values():
                line.set_data([], [])
            for line in est_lines.values():
                line.set_data([], [])
            meas_scatter.set_offsets(np.empty((0, 2)))
            return list(true_lines.values()) + list(est_lines.values()) + [meas_scatter]
        
        def update(frame):
            # 更新真实轨迹
            for target_id, trajectory in true_trajectories.items():
                # 找到该时刻之前的数据
                valid_idx = frame + 1
                if valid_idx > len(trajectory):
                    valid_idx = len(trajectory)
                true_lines[target_id].set_data(trajectory[:valid_idx, 0], 
                                               trajectory[:valid_idx, 1])
            
            # 更新估计轨迹
            if estimated_trajectories:
                for target_id, trajectory in estimated_trajectories.items():
                    if target_id in est_lines:
                        valid_idx = min(frame + 1, len(trajectory))
                        est_lines[target_id].set_data(trajectory[:valid_idx, 0],
                                                     trajectory[:valid_idx, 1])
            
            # 更新观测
            if frame < len(measurements):
                meas = measurements[frame]
                if meas is not None and len(meas) > 0:
                    meas_scatter.set_offsets(meas)
                else:
                    meas_scatter.set_offsets(np.empty((0, 2)))
            
            ax.set_title(f'Tracking Animation - Step {frame}')
            
            return list(true_lines.values()) + list(est_lines.values()) + [meas_scatter]
        
        anim = animation.FuncAnimation(fig, update, frames=n_steps,
                                       init_func=init, blit=True, interval=interval)
        
        if save_path:
            anim.save(save_path, writer='pillow', fps=1000//interval)
        
        return anim
    
    def plot_target_positions(self,
                              true_positions: Dict[int, np.ndarray],
                              estimated_positions: Dict[int, np.ndarray],
                              time_steps: np.ndarray,
                              target_id: int,
                              title: Optional[str] = None,
                              save_path: Optional[str] = None) -> plt.Figure:
        """绘制单个目标的位置对比
        
        Args:
            true_positions: 真实位置字典
            estimated_positions: 估计位置字典
            time_steps: 时间步数组
            target_id: 目标ID
            title: 图形标题
            save_path: 保存路径
            
        Returns:
            matplotlib图形对象
        """
        fig, axes = plt.subplots(2, 1, figsize=self.figsize, sharex=True)
        
        true_pos = true_positions.get(target_id)
        est_pos = estimated_positions.get(target_id)
        
        if true_pos is not None:
            axes[0].plot(time_steps[:len(true_pos)], true_pos[:, 0], 
                        'b-', linewidth=2, label='True')
            axes[1].plot(time_steps[:len(true_pos)], true_pos[:, 1], 
                        'b-', linewidth=2, label='True')
        
        if est_pos is not None:
            axes[0].plot(time_steps[:len(est_pos)], est_pos[:, 0], 
                        'r--', linewidth=1.5, label='Estimated')
            axes[1].plot(time_steps[:len(est_pos)], est_pos[:, 1], 
                        'r--', linewidth=1.5, label='Estimated')
        
        axes[0].set_ylabel('X Position')
        axes[1].set_ylabel('Y Position')
        axes[1].set_xlabel('Time')
        
        if title:
            fig.suptitle(title)
        else:
            fig.suptitle(f'Target {target_id} Position')
        
        for ax in axes:
            ax.legend()
            ax.grid(True, alpha=0.3)
        
        plt.tight_layout()
        
        if save_path:
            plt.savefig(save_path, dpi=300, bbox_inches='tight')
        
        return fig
    
    def plot_ospa_over_time(self,
                            ospa_values: np.ndarray,
                            localization_errors: np.ndarray,
                            cardinality_errors: np.ndarray,
                            time_steps: np.ndarray,
                            title: str = "OSPA Error Over Time",
                            save_path: Optional[str] = None) -> plt.Figure:
        """绘制OSPA误差随时间变化（参考MATLAB RFS Toolbox）
        
        上图：OSPA总误差
        下图：定位误差 + 基数误差分解
        
        Args:
            ospa_values: OSPA值数组
            localization_errors: 定位误差数组
            cardinality_errors: 基数误差数组
            time_steps: 时间步数组
            title: 图形标题
            save_path: 保存路径
            
        Returns:
            matplotlib图形对象
        """
        fig, axes = plt.subplots(2, 1, figsize=(12, 8), sharex=True)
        
        # 上图：OSPA总误差
        axes[0].plot(time_steps, ospa_values, 'b-', linewidth=1.5, label='OSPA')
        axes[0].set_ylabel('OSPA Distance')
        axes[0].set_title(title)
        axes[0].legend(loc='best')
        axes[0].grid(True, alpha=0.3)
        
        # 下图：误差分解
        axes[1].plot(time_steps, localization_errors, 'g-', linewidth=1.5, label='Localization')
        axes[1].plot(time_steps, cardinality_errors, 'r-', linewidth=1.5, label='Cardinality')
        axes[1].set_xlabel('Time Step')
        axes[1].set_ylabel('Error Component')
        axes[1].legend(loc='best')
        axes[1].grid(True, alpha=0.3)
        
        plt.tight_layout()
        
        if save_path:
            plt.savefig(save_path, dpi=300, bbox_inches='tight')
        
        return fig
    
    def plot_cardinality_over_time(self,
                                   true_cardinality: np.ndarray,
                                   estimated_cardinality: np.ndarray,
                                   time_steps: np.ndarray,
                                   title: str = "Cardinality Estimation",
                                   save_path: Optional[str] = None) -> plt.Figure:
        """绘制基数估计随时间变化
        
        Args:
            true_cardinality: 真实基数数组
            estimated_cardinality: 估计基数数组
            time_steps: 时间步数组
            title: 图形标题
            save_path: 保存路径
            
        Returns:
            matplotlib图形对象
        """
        fig, ax = plt.subplots(figsize=(12, 4))
        
        ax.plot(time_steps, true_cardinality, 'b-', linewidth=2, label='True')
        ax.plot(time_steps, estimated_cardinality, 'r--', linewidth=1.5, label='Estimated')
        ax.set_xlabel('Time Step')
        ax.set_ylabel('Number of Targets')
        ax.set_title(title)
        ax.legend(loc='best')
        ax.grid(True, alpha=0.3)
        ax.set_ylim(bottom=-0.5)
        
        plt.tight_layout()
        
        if save_path:
            plt.savefig(save_path, dpi=300, bbox_inches='tight')
        
        return fig
