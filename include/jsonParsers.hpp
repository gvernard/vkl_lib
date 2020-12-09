#ifndef JSON_PARSERS_HPP
#define JSON_PARSERS_HPP

#include "json/json.h"

class CollectionMassModels;
class CollectionProfiles;

class JsonParsers {
public:
  static CollectionMassModels parse_mass_model(const Json::Value mass_model);
  static CollectionProfiles parse_profile(const Json::Value profile);
};

#endif /* JSON_PARSERS_HPP */
