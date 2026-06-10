"""
贝叶斯滤波示例项目 - 主运行脚本
"""
import sys
import os
from pathlib import Path

# 添加项目根目录到Python路径
PROJECT_ROOT = Path(__file__).parent
sys.path.insert(0, str(PROJECT_ROOT))

import argparse

from demo import (
    run_linear_tracking_demo,
    run_nonlinear_tracking_demo,
    run_multi_target_demo,
    run_performance_comparison_demo
)


def main():
    """主函数"""
    parser = argparse.ArgumentParser(description='贝叶斯滤波示例项目')
    
    parser.add_argument('--demo', type=str, default='linear',
                       choices=['linear', 'nonlinear', 'multi', 'compare'],
                       help='演示类型: linear, nonlinear, multi, compare')
    
    parser.add_argument('--duration', type=float, default=50.0,
                       help='场景时长（秒）')
    
    parser.add_argument('--n_targets', type=int, default=1,
                       help='目标数量')
    
    parser.add_argument('--filter', type=str, default='UKF',
                       choices=['EKF', 'UKF', 'CKF'],
                       help='非线性滤波器类型')
    
    parser.add_argument('--association', type=str, default='knn',
                       choices=['nn', 'knn', 'jpda', 'mht', 'phd'],
                       help='数据关联算法 (仅 multi 模式): nn, knn, jpda, mht, phd')
    
    parser.add_argument('--no_plot', action='store_true',
                       help='不显示图形')
    
    parser.add_argument('--save', type=str, default=None,
                       help='保存路径')
    
    args = parser.parse_args()
    
    print("=" * 60)
    print("贝叶斯滤波示例项目")
    print("=" * 60)
    
    if args.demo == 'linear':
        print("\n运行线性跟踪演示...")
        results = run_linear_tracking_demo(
            duration=args.duration,
            n_targets=args.n_targets,
            show_plot=not args.no_plot,
            save_path=args.save
        )
        
    elif args.demo == 'nonlinear':
        print(f"\n运行非线性跟踪演示 ({args.filter})...")
        results = run_nonlinear_tracking_demo(
            duration=args.duration,
            n_targets=args.n_targets,
            filter_type=args.filter,
            show_plot=not args.no_plot,
            save_path=args.save
        )
        
    elif args.demo == 'multi':
        print(f"\n运行多目标跟踪演示 (关联算法: {args.association.upper()})...")
        results = run_multi_target_demo(
            duration=args.duration,
            n_targets=args.n_targets,
            association_method=args.association,
            show_plot=not args.no_plot,
            save_path=args.save
        )
        
    elif args.demo == 'compare':
        print("\n运行性能对比演示...")
        results = run_performance_comparison_demo(
            duration=args.duration,
            n_targets=args.n_targets,
            show_plot=not args.no_plot,
            save_path=args.save
        )
    
    else:
        print(f"未知的演示类型: {args.demo}")
        return
    
    print("\n" + "=" * 60)
    print("演示完成!")
    print("=" * 60)


if __name__ == '__main__':
    main()
