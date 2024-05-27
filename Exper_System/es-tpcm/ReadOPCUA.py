from opcua import Client

def read_opcua_values(server_url, node_ids):
    # Create a client
    client = Client(server_url)

    try:
        # Connect to the OPC UA server
        client.connect()

        # Create a list to store the read values
        values = []

        # Read values for each specified NodeID
        for node_id in node_ids:
            # Find the Node object by NodeID
            node = client.get_node(node_id)

            # Read the value attribute of the Node
            value = node.get_value()
            values.append((node_id, value))

        return values

    finally:
        # Disconnect from the OPC UA server
        client.disconnect()

if __name__ == "__main__":
    # Specify the OPC UA server URL
    server_url = "opc.tcp://000.000.000.000:4840"

    # Specify the NodeIDs you want to read
    node_ids = [
        "ns=4;s=;s=MAIN.VarPOU_LoTuS.LoTuS.DAK.BvL_Tanks.BvL_Tank_RZ.Tankheizung.BvL_T_T.sensorState.fValue",
    ]

    # Read OPC UA values
    result = read_opcua_values(server_url, node_ids)

    # Display the results
    for node_id, value in result:
        print(f"NodeID: {node_id}, Value: {value}")
