import pytest
import json
import tempfile
import os
import pandas as pd
import numpy as np
from aiohttp import web
from aiohttp.test_utils import TestClient, setup_test_loop

from main import (
    BedTileRenderer,
    BigWigRenderer,
    HiCRangeRenderer,
    GlobalManager,
    main,
    handle_register,
    handle_health,
    bed_info,
    bed_tiles,
    handle_1d_tile,
    handle_mcool_chroms,
    hic_info,
    hic_range,
    handle_hic_tile,
    bed_global_info,
    bed_get_prediction,
    handle_tad_status,
    handle_get_all_records,
    handle_delete_record,
    handle_delete_all_records,
    get_genome_index
)


@pytest.fixture
def test_bed_file():
    """创建测试用的BED文件"""
    with tempfile.NamedTemporaryFile(mode='w', suffix='.bed', delete=False) as f:
        f.write("""1\t10000\t20000\tregion1
1\t30000\t40000\tregion2
2\t50000\t60000\tregion3
2\t70000\t80000\tregion4
""")
    yield f.name
    os.unlink(f.name)


@pytest.fixture
def test_app():
    """创建测试用的aiohttp应用"""
    app = web.Application(client_max_size=1024 ** 2 * 100)
    app['manager'] = GlobalManager()
    
    app.add_routes([
        web.get("/api/health", handle_health),
        web.post("/api/register", handle_register),
        web.get('/api/genome/index/{token}', get_genome_index),
        web.get("/api/records", handle_get_all_records),
        web.get('/api/bed/info/{token}/{chrom}', bed_info),
        web.get('/api/bed/tiles/{token}/{chrom}/{zoom}/{x}', bed_tiles),
        web.delete("/api/records/{token}", handle_delete_record),
        web.delete("/api/records", handle_delete_all_records),
        web.get('/api/v1/tiles/1d/{token}/{z}/{x}', handle_1d_tile),
        web.get('/api/bed/global-info/{token}', bed_global_info),
        web.get('/api/bed/prediction/{token}', bed_get_prediction),
        web.get("/api/mcool/chroms/{token}", handle_mcool_chroms),
        web.get("/api/hic/info/{token}/{chrom}", hic_info),
        web.get("/api/hic/range/{token}/{chrom}", hic_range),
        web.get('/api/v1/tiles/2d/{token}/{z}/{x}/{y}', handle_hic_tile),
        web.get("/api/tad/status/{token}", handle_tad_status),
    ])
    
    return app


@pytest.mark.asyncio
async def test_bed_renderer_loading(test_bed_file):
    """测试BED渲染器正确加载和解析文件"""
    renderer = BedTileRenderer(test_bed_file)
    
    assert len(renderer.df) == 4
    assert set(renderer.chrom_groups.keys()) == {'1', '2'}
    
    # 测试染色体归一化
    assert renderer._normalize_chrom('chr1') == '1'
    assert renderer._normalize_chrom('Chr2') == '2'
    assert renderer._normalize_chrom('3') == '3'
    
    # 测试染色体信息
    chroms = renderer.get_all_chroms()
    assert len(chroms) == 2
    assert chroms[0]['name'] == '1'
    assert chroms[0]['length'] == 40000
    assert chroms[1]['name'] == '2'
    assert chroms[1]['length'] == 80000


@pytest.mark.asyncio
async def test_bed_renderer_fetch_tile(test_bed_file):
    """测试BED瓦片获取功能"""
    renderer = BedTileRenderer(test_bed_file)
    
    # 测试正常缩放级别
    data = await renderer.fetch_tile('1', 0, 0)
    assert len(data) == 2
    
    # 测试不存在的染色体
    data = await renderer.fetch_tile('99', 0, 0)
    assert data == []
    
    # 测试超出缩放范围
    data = await renderer.fetch_tile('1', 100, 0)
    assert len(data) > 0


@pytest.mark.asyncio
async def test_bed_info_endpoint(test_app, test_bed_file):
    """测试BED信息API端点"""
    async with TestClient(test_app) as client:
        # 先注册文件
        resp = await client.post('/api/register', json={
            'path': test_bed_file,
            'type': 'bed'
        })
        assert resp.status == 200
        token = (await resp.json())['token']
        
        # 获取染色体信息
        resp = await client.get(f'/api/bed/info/{token}/1')
        assert resp.status == 200
        data = await resp.json()
        assert data['chromosome'] == '1'
        assert 'available_resolutions' in data
        assert 'max_zoom' in data
        
        # 测试不存在的染色体
        resp = await client.get(f'/api/bed/info/{token}/99')
        assert resp.status == 200
        assert 'error' in await resp.json()


@pytest.mark.asyncio
async def test_bed_tiles_endpoint(test_app, test_bed_file):
    """测试BED瓦片API端点"""
    async with TestClient(test_app) as client:
        resp = await client.post('/api/register', json={
            'path': test_bed_file,
            'type': 'bed'
        })
        token = (await resp.json())['token']
        
        resp = await client.get(f'/api/bed/tiles/{token}/1/0/0')
        assert resp.status == 200
        data = await resp.json()
        assert 'data' in data
        assert 'zoom' in data
        assert 'x' in data


@pytest.mark.asyncio
async def test_1d_tile_endpoint(test_app, test_bed_file):
    """测试1D瓦片API端点"""
    async with TestClient(test_app) as client:
        resp = await client.post('/api/register', json={
            'path': test_bed_file,
            'type': 'bed'
        })
        token = (await resp.json())['token']
        
        resp = await client.get(f'/api/v1/tiles/1d/{token}/0/0')
        assert resp.status == 200
        data = await resp.json()
        assert 'tile_id' in data
        assert 'type' in data
        assert data['type'] == 'bed'
        assert 'data' in data


@pytest.mark.asyncio
async def test_bed_global_info(test_app, test_bed_file):
    """测试全局BED信息端点"""
    async with TestClient(test_app) as client:
        resp = await client.post('/api/register', json={
            'path': test_bed_file,
            'type': 'bed'
        })
        token = (await resp.json())['token']
        
        resp = await client.get(f'/api/bed/global-info/{token}')
        assert resp.status == 200
        data = await resp.json()
        assert 'genome_length' in data
        assert 'available_resolutions' in data
        assert 'chromosomes' in data
        assert len(data['chromosomes']) == 2


@pytest.mark.asyncio
async def test_health_endpoint(test_app):
    """测试健康检查端点"""
    async with TestClient(test_app) as client:
        resp = await client.get('/api/health')
        assert resp.status == 200
        data = await resp.json()
        assert data['status'] == 'ok'
        assert data['service'] == 'genomics-api'


@pytest.mark.asyncio
async def test_register_endpoint(test_app, test_bed_file):
    """测试文件注册端点"""
    async with TestClient(test_app) as client:
        # 正常注册
        resp = await client.post('/api/register', json={
            'path': test_bed_file,
            'type': 'bed'
        })
        assert resp.status == 200
        data = await resp.json()
        assert 'token' in data
        assert data['type'] == 'bed'
        
        # 缺少path参数
        resp = await client.post('/api/register', json={})
        assert resp.status == 400
        
        # 无效JSON
        resp = await client.post('/api/register', data='invalid json')
        assert resp.status == 400


@pytest.mark.asyncio
async def test_genome_index_endpoint(test_app, test_bed_file):
    """测试基因组索引端点"""
    async with TestClient(test_app) as client:
        resp = await client.post('/api/register', json={
            'path': test_bed_file,
            'type': 'bed'
        })
        token = (await resp.json())['token']
        
        resp = await client.get(f'/api/genome/index/{token}')
        assert resp.status == 200
        data = await resp.json()
        assert 'total_length' in data
        assert 'chromosomes' in data


@pytest.mark.asyncio
async def test_global_manager():
    """测试全局管理器功能"""
    manager = GlobalManager()
    
    # 测试空查询
    assert manager.get('invalid-token') is None
    
    # 测试注册重复文件返回相同token
    with tempfile.NamedTemporaryFile(mode='w', suffix='.bed', delete=False) as f:
        f.write("1\t1000\t2000\ttest\n")
    
    token1 = manager.register(f.name, 'bed')
    token2 = manager.register(f.name, 'bed')
    assert token1 == token2
    
    # 测试获取绝对信息
    offsets, total_len, chrom_list = manager.get_abs_info(token1)
    assert isinstance(offsets, dict)
    assert total_len > 0
    assert len(chrom_list) > 0
    
    # 测试坐标转换
    chrom, pos = manager.abs_to_rel(1500, offsets)
    assert chrom == '1'
    assert pos == 1500
    
    # 关闭所有资源
    manager.close_all()
    
    os.unlink(f.name)


def test_bed_renderer_edge_cases():
    """测试BED渲染器边缘情况"""
    # 测试只有3列的BED文件
    with tempfile.NamedTemporaryFile(mode='w', suffix='.bed', delete=False) as f:
        f.write("""1\t10000\t20000
1\t30000\t40000
""")
    
    renderer = BedTileRenderer(f.name)
    assert len(renderer.df) == 2
    assert 'name' in renderer.df.columns
    assert all(renderer.df['name'] == "")
    
    os.unlink(f.name)


def test_bed_renderer_prediction_loading():
    """测试预测概率加载功能"""
    # 创建BED文件
    with tempfile.NamedTemporaryFile(mode='w', suffix='.bed', delete=False) as f:
        f.write("1\t10000\t20000\tregion1\n")
    bed_path = f.name
    
    # 创建预测文件
    with tempfile.NamedTemporaryFile(mode='w', suffix='.txt', delete=False) as f:
        f.write("chrom\tstart\tend\tprediction\n")
        f.write("1\t10000\t20000\t0.85\n")
    pred_path = f.name
    
    renderer = BedTileRenderer(bed_path, prediction_path=pred_path)
    assert len(renderer.prediction_map) == 1
    assert renderer.prediction_map[('1', 10000, 20000)] == 0.85
    
    os.unlink(bed_path)
    os.unlink(pred_path)


@pytest.mark.asyncio
async def test_bed_prediction_endpoint(test_app):
    """测试BED预测端点"""
    # 创建带预测数据的BED渲染器
    with tempfile.NamedTemporaryFile(mode='w', suffix='.bed', delete=False) as f:
        f.write("1\t10000\t20000\tregion1\n")
    bed_path = f.name
    
    with tempfile.NamedTemporaryFile(mode='w', suffix='.txt', delete=False) as f:
        f.write("chrom\tstart\tend\tprediction\n")
        f.write("1\t10000\t20000\t0.85\n")
    pred_path = f.name
    
    async with TestClient(test_app) as client:
        resp = await client.post('/api/register', json={
            'path': bed_path,
            'type': 'bed',
            'prediction_path': pred_path
        })
        token = (await resp.json())['token']
        
        # 测试存在的区间
        resp = await client.get(f'/api/bed/prediction/{token}?chrom=1&start=10000&end=20000')
        assert resp.status == 200
        data = await resp.json()
        assert data['prediction'] == 0.85
        
        # 测试不存在的区间
        resp = await client.get(f'/api/bed/prediction/{token}?chrom=1&start=99999&end=999999')
        assert resp.status == 200
        data = await resp.json()
        assert data['prediction'] is None
        
        # 测试无效坐标
        resp = await client.get(f'/api/bed/prediction/{token}?chrom=1&start=abc&end=def')
        assert resp.status == 400
    
    os.unlink(bed_path)
    os.unlink(pred_path)


def test_invalid_bed_file():
    """测试无效BED文件处理"""
    with tempfile.NamedTemporaryFile(mode='w', suffix='.bed', delete=False) as f:
        f.write("invalid content\n")
    
    with pytest.raises(Exception):
        BedTileRenderer(f.name)
    
    os.unlink(f.name)


if __name__ == "__main__":
    pytest.main([__file__, "-v"])