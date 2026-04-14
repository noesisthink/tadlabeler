import asyncio
import aiohttp
import json

BASE_URL = "http://127.0.0.1:5001"


BIGWIG_PATH = "test.bw"


async def register_bigwig(session):
    """注册 BigWig 文件"""
    url = f"{BASE_URL}/api/register"

    payload = {
        "path": BIGWIG_PATH,
        "type": "bigwig"
    }

    async with session.post(url, json=payload) as resp:
        data = await resp.json()
        print("\n[REGISTER RESPONSE]")
        print(json.dumps(data, indent=2))

        return data.get("token")


async def get_genome_index(session, token):
    """获取全基因组索引"""
    url = f"{BASE_URL}/api/genome/index/{token}"

    async with session.get(url) as resp:
        data = await resp.json()
        print("\n[GENOME INDEX]")
        print(json.dumps(data, indent=2))


async def fetch_tile(session, token, z=0, x=0):
    """获取 BigWig tile"""
    url = f"{BASE_URL}/api/v1/tiles/1d/{token}/{z}/{x}"

    async with session.get(url) as resp:
        data = await resp.json()

        print(f"\n[TILE z={z}, x={x}]")
        print("type:", data.get("type"))
        print("abs_range:", data.get("abs_range"))
        print("data (前20个):", data.get("data")[:20])


async def main():
    async with aiohttp.ClientSession() as session:
        # 1. 注册
        token = await register_bigwig(session)
        if not token:
            print("❌ 注册失败")
            return

        # 2. genome index
        await get_genome_index(session, token)

        for z in range(0, 3):
            for x in range(0, 2):
                await fetch_tile(session, token, z, x)


if __name__ == "__main__":
    asyncio.run(main())