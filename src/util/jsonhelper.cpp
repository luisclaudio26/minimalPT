#include "../../include/util/jsonhelper.h"
#include <glm/gtx/transform.hpp>

namespace JSONHelper
{
	Vec4 load_vec4(const nlohmann::json& in)
	{
		Vec4 out;
		out.x = in[0].get<double>();
		out.y = in[1].get<double>();
		out.z = in[2].get<double>();
		out.w = in[3].get<double>();
		return out;
	}

	Vec3 load_vec3(const nlohmann::json& in)
	{
		Vec3 out;
		out.x = in[0].get<double>();
		out.y = in[1].get<double>();
		out.z = in[2].get<double>();
		return out;
	}

	Vec2 load_vec2(const nlohmann::json& in)
	{
		Vec2 out;
		out.x = in[0].get<double>();
		out.y = in[1].get<double>();
		return out;
	}

	Mat4 load_transformation(const nlohmann::json& in)
	{
		Mat4 out(1.0f);
		Mat4 scl(1.0f), rot(1.0f), trans(1.0f);

		if( in.find("scale") != in.end() )
		{
			glm::vec3 scl_v = JSONHelper::load_vec3( in["scale"] );
			scl = glm::scale(scl_v);
		}

		if( in.find("rotate") != in.end() )
		{
			glm::vec4 rot_v = JSONHelper::load_vec4( in["rotate"] );
			rot = glm::rotate(rot_v[3], glm::vec3(rot_v));
		}

		if( in.find("translate") != in.end() )
		{
			glm::vec3 trans_v = JSONHelper::load_vec3( in["translate"] );
			trans = glm::translate(trans_v);
		}

		std::string order = in["order"].get<std::string>();
		for(auto t = order.begin(); t != order.end(); ++t)
		{
			if(*t == 't' || *t == 'T')
				out = trans * out;
			if(*t == 'r' || *t == 'R')
				out = rot * out;
			if(*t == 's' || *t == 'S')
				out = scl * out;
		}

		return out;
	}

}
